/*

   Native PDBx/mmCIF reader for the rigid receptor.

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.

   The receptor atoms are read straight from the _atom_site category; no intermediate PDBQT is
   produced. PDBQT carries AutoDock atom types that mmCIF does not, so they are derived here:

   - element (type_symbol) gives the base type; C/N/O/S/H are refined below.
   - covalent connectivity is perceived with the same criterion as model::assign_bonds.
   - residues that contain hydrogens are typed from connectivity alone, like a prepared PDBQT:
     H bonded to N/O/S -> HD, N without H and with <= 2 non-metal neighbours -> NA (not for
     backbone N), O -> OA, S -> SA.
     Donors then come out of model::assign_types through the HD atoms, exactly as for PDBQT.
   - residues without hydrogens (typical for PDB/AlphaFold files) use built-in templates of
     standard amino acids, nucleotides and water to know which heavy atoms carry polar
     hydrogens; histidine tautomers come from their hydrogen-bond partners, and residues
     without a template (ligands, cofactors, modified residues) are typed from bond lengths,
     angles and ring planarity. Those atoms are reported in rigid::implicit_donors and promoted
     to XS donors after typing (see parse_receptor_mmcif).

*/

#include "parse_mmcif.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdlib>
#include <initializer_list>
#include <iostream>
#include <map>
#include <memory>
#include <set>
#include <sstream>
#include <unordered_map>
#include "atom_constants.h"
#include "file.h"
#include "parse_error.h"

namespace {

std::string to_lower(std::string s) {
    for (char& ch : s) ch = char(std::tolower(static_cast<unsigned char>(ch)));
    return s;
}

std::string to_upper(std::string s) {
    for (char& ch : s) ch = char(std::toupper(static_cast<unsigned char>(ch)));
    return s;
}

bool ends_with_ci(const std::string& str, const std::string& suffix) {
    return str.size() >= suffix.size()
           && to_lower(str.substr(str.size() - suffix.size())) == suffix;
}

// ------------------------------------------------------------------------------------------------
// CIF tokenizer (STAR/CIF 1.1 lexical rules: bare words, '...' and "..." strings that close only
// before whitespace, ;-delimited text fields, # comments)

struct cif_token {
    std::string text;
    bool quoted = false;  // quoted strings and text fields are never keywords, tags or nulls
    unsigned line = 0;

    bool is_tag() const { return !quoted && !text.empty() && text[0] == '_'; }
    bool is_keyword() const {
        if (quoted) return false;
        const std::string l = to_lower(text);
        return starts_with(l, "data_") || starts_with(l, "save_") || l == "loop_"
               || l == "global_" || l == "stop_";
    }
    bool is_null() const { return !quoted && (text == "." || text == "?"); }
};

class cif_tokenizer {
public:
    explicit cif_tokenizer(std::istream& in) : m_in(in) {}

    bool next(cif_token& tok) {
        if (m_has_pushed) {
            tok = m_pushed;
            m_has_pushed = false;
            return true;
        }
        while (true) {
            if (m_pos >= m_line.size()) {
                if (!fetch_line()) return false;
                if (!m_line.empty() && m_line[0] == ';') {
                    read_text_field(tok);
                    return true;
                }
            }
            while (m_pos < m_line.size() && std::isspace(static_cast<unsigned char>(m_line[m_pos])))
                ++m_pos;
            if (m_pos >= m_line.size()) continue;

            const char ch = m_line[m_pos];
            if (ch == '#') {
                m_pos = m_line.size();
                continue;
            }
            tok.line = m_line_no;
            if (ch == '\'' || ch == '"') {
                sz end = m_pos + 1;
                while (end < m_line.size()
                       && !(m_line[end] == ch
                            && (end + 1 == m_line.size()
                                || std::isspace(static_cast<unsigned char>(m_line[end + 1])))))
                    ++end;
                if (end >= m_line.size())
                    throw struct_parse_error(
                        "Unterminated quoted string in mmCIF file (line " + std::to_string(m_line_no)
                            + ").",
                        m_line);
                tok.text = m_line.substr(m_pos + 1, end - m_pos - 1);
                tok.quoted = true;
                m_pos = end + 1;
                return true;
            }
            sz end = m_pos;
            while (end < m_line.size() && !std::isspace(static_cast<unsigned char>(m_line[end])))
                ++end;
            tok.text = m_line.substr(m_pos, end - m_pos);
            tok.quoted = false;
            m_pos = end;
            return true;
        }
    }

    void push_back(const cif_token& tok) {
        m_pushed = tok;
        m_has_pushed = true;
    }

private:
    bool fetch_line() {
        m_pos = 0;
        if (!std::getline(m_in, m_line)) {
            m_line.clear();
            return false;
        }
        if (!m_line.empty() && m_line.back() == '\r') m_line.pop_back();
        ++m_line_no;
        return true;
    }

    void read_text_field(cif_token& tok) {
        tok.line = m_line_no;
        tok.quoted = true;
        tok.text = m_line.substr(1);
        while (true) {
            if (!fetch_line())
                throw struct_parse_error("Unterminated text field (;) starting at line "
                                         + std::to_string(tok.line) + " of mmCIF file.");
            if (!m_line.empty() && m_line[0] == ';') {
                m_pos = 1;
                return;
            }
            tok.text += '\n';
            tok.text += m_line;
        }
    }

    std::istream& m_in;
    std::string m_line;
    sz m_pos = 0;
    unsigned m_line_no = 0;
    cif_token m_pushed;
    bool m_has_pushed = false;
};

// ------------------------------------------------------------------------------------------------
// _atom_site rows

struct atom_site_columns {
    int type_symbol = -1;
    int label_atom_id = -1, auth_atom_id = -1;
    int label_comp_id = -1, auth_comp_id = -1;
    int label_asym_id = -1, auth_asym_id = -1;
    int label_seq_id = -1, auth_seq_id = -1;
    int ins_code = -1, alt_id = -1, model_num = -1;
    int x = -1, y = -1, z = -1;

    explicit atom_site_columns(const std::vector<std::string>& tags) {
        VINA_FOR_IN(i, tags) {
            const std::string& t = tags[i];
            const int idx = int(i);
            if (t == "_atom_site.type_symbol") type_symbol = idx;
            else if (t == "_atom_site.label_atom_id") label_atom_id = idx;
            else if (t == "_atom_site.auth_atom_id") auth_atom_id = idx;
            else if (t == "_atom_site.label_comp_id") label_comp_id = idx;
            else if (t == "_atom_site.auth_comp_id") auth_comp_id = idx;
            else if (t == "_atom_site.label_asym_id") label_asym_id = idx;
            else if (t == "_atom_site.auth_asym_id") auth_asym_id = idx;
            else if (t == "_atom_site.label_seq_id") label_seq_id = idx;
            else if (t == "_atom_site.auth_seq_id") auth_seq_id = idx;
            else if (t == "_atom_site.pdbx_pdb_ins_code") ins_code = idx;
            else if (t == "_atom_site.label_alt_id") alt_id = idx;
            else if (t == "_atom_site.pdbx_pdb_model_num") model_num = idx;
            else if (t == "_atom_site.cartn_x") x = idx;
            else if (t == "_atom_site.cartn_y") y = idx;
            else if (t == "_atom_site.cartn_z") z = idx;
        }
        if (x < 0 || y < 0 || z < 0)
            throw struct_parse_error(
                "mmCIF _atom_site is missing Cartn_x, Cartn_y or Cartn_z.");
        if (label_atom_id < 0 && auth_atom_id < 0)
            throw struct_parse_error("mmCIF _atom_site is missing label_atom_id/auth_atom_id.");
    }
};

// nulls ('.' and '?') are stored as empty strings
typedef std::vector<std::string> cif_row;

const std::string& field(const cif_row& row, int idx) {
    static const std::string empty;
    return idx < 0 ? empty : row[sz(idx)];
}

const std::string& first_of(const cif_row& row, int preferred, int fallback) {
    const std::string& v = field(row, preferred);
    return v.empty() ? field(row, fallback) : v;
}

struct raw_atom {
    std::string element;  // capitalised, e.g. "C", "Zn"
    std::string name;     // atom name, '*' in old nucleotide names converted to '\''
    std::string res_name;
    std::string res_key;
    vec coords;
    unsigned line;
};

fl parse_coordinate(const std::string& s, unsigned line) {
    const char* begin = s.c_str();
    char* end = nullptr;
    const double v = std::strtod(begin, &end);
    // allow a trailing standard uncertainty, e.g. "12.345(6)"
    if (s.empty() || end == begin || (*end != '\0' && *end != '(') || !std::isfinite(v))
        throw struct_parse_error("Coordinate \"" + s + "\" is not valid (mmCIF line "
                                 + std::to_string(line) + ").");
    return fl(v);
}

std::string normalize_element(const std::string& type_symbol, const std::string& atom_name,
                              const std::string& res_name) {
    std::string letters;
    if (!type_symbol.empty()) {
        for (char ch : type_symbol) {  // "FE2+" -> "FE"
            if (!std::isalpha(static_cast<unsigned char>(ch))) break;
            letters += ch;
        }
    } else {  // no type_symbol column: PDB-style guess from the atom name
        std::string name;
        for (char ch : atom_name)
            if (std::isalpha(static_cast<unsigned char>(ch))) name += ch;
        if (name.size() >= 2 && to_upper(name) == to_upper(res_name))
            letters = name.substr(0, 2);  // single-atom ions: ZN/ZN, MG/MG
        else
            letters = name.substr(0, 1);
    }
    if (letters.empty()) return letters;
    std::string el(1, char(std::toupper(static_cast<unsigned char>(letters[0]))));
    el += to_lower(letters.substr(1));
    return el;
}

class atom_site_reader {
public:
    std::vector<raw_atom> atoms;
    std::set<std::string> skipped_models;
    sz skipped_altloc_atoms = 0;

    void add_row(const atom_site_columns& col, const cif_row& row, unsigned line) {
        const std::string& model = field(row, col.model_num);
        if (!m_have_model) {
            m_model = model;
            m_have_model = true;
        }
        if (model != m_model) {
            skipped_models.insert(model);
            return;
        }

        raw_atom a;
        a.name = first_of(row, col.label_atom_id, col.auth_atom_id);
        std::replace(a.name.begin(), a.name.end(), '*', '\'');
        a.res_name = first_of(row, col.auth_comp_id, col.label_comp_id);
        const std::string& chain = first_of(row, col.label_asym_id, col.auth_asym_id);
        const std::string& seq = first_of(row, col.auth_seq_id, col.label_seq_id);
        const std::string& ins = field(row, col.ins_code);
        a.res_key = chain + '\x1f' + seq + '\x1f' + ins + '\x1f' + a.res_name;

        // alternate locations: keep the first one met in each residue
        const std::string& alt = field(row, col.alt_id);
        if (!alt.empty()) {
            std::map<std::string, std::string>::iterator it = m_residue_alt.find(a.res_key);
            if (it == m_residue_alt.end())
                m_residue_alt[a.res_key] = alt;
            else if (it->second != alt) {
                ++skipped_altloc_atoms;
                return;
            }
        }

        a.element = normalize_element(field(row, col.type_symbol), a.name, a.res_name);
        a.coords = vec(parse_coordinate(field(row, col.x), line),
                       parse_coordinate(field(row, col.y), line),
                       parse_coordinate(field(row, col.z), line));
        a.line = line;
        atoms.push_back(a);
    }

private:
    bool m_have_model = false;
    std::string m_model;
    std::map<std::string, std::string> m_residue_alt;
};

void read_atom_site(std::istream& in, atom_site_reader& reader) {
    cif_tokenizer tz(in);
    cif_token tok;
    bool in_block = false;
    bool found = false;
    std::vector<std::string> single_tags;  // non-looped _atom_site (a single atom)
    cif_row single_row;
    unsigned single_line = 0;

    while (tz.next(tok)) {
        const std::string lower = to_lower(tok.text);
        if (!tok.quoted && starts_with(lower, "data_")) {
            if (in_block) break;  // only the first data block describes the receptor
            in_block = true;
        } else if (!tok.quoted && lower == "loop_") {
            std::vector<std::string> tags;
            bool more;
            while ((more = tz.next(tok)) && tok.is_tag()) tags.push_back(to_lower(tok.text));
            const bool is_atom_site
                = !tags.empty() && starts_with(tags.front(), "_atom_site.");
            std::unique_ptr<atom_site_columns> col;
            if (is_atom_site) col.reset(new atom_site_columns(tags));

            cif_row row;
            unsigned row_line = 0;
            while (more && !tok.is_tag() && !tok.is_keyword()) {
                if (is_atom_site) {
                    if (row.empty()) row_line = tok.line;
                    row.push_back(tok.is_null() ? std::string() : tok.text);
                    if (row.size() == tags.size()) {
                        reader.add_row(*col, row, row_line);
                        row.clear();
                    }
                }
                more = tz.next(tok);
            }
            if (is_atom_site) {
                if (!row.empty())
                    throw struct_parse_error(
                        "mmCIF _atom_site loop has a number of values that is not a multiple of "
                        "its number of columns (line " + std::to_string(row_line) + ").");
                found = true;
                break;
            }
            if (more) tz.push_back(tok);
        } else if (tok.is_tag()) {
            const std::string tag = lower;
            if (!tz.next(tok) || tok.is_tag() || tok.is_keyword())
                throw struct_parse_error("mmCIF tag " + tag + " has no value (line "
                                         + std::to_string(tok.line) + ").");
            if (starts_with(tag, "_atom_site.")) {
                if (single_tags.empty()) single_line = tok.line;
                single_tags.push_back(tag);
                single_row.push_back(tok.is_null() ? std::string() : tok.text);
            }
        }
        // anything else (save frames, global_, stray values) is irrelevant for _atom_site
    }

    if (!found && !single_tags.empty()) {
        atom_site_columns col(single_tags);
        reader.add_row(col, single_row, single_line);
        found = true;
    }
    if (!found) throw struct_parse_error("No _atom_site category found in mmCIF file.");
    if (reader.atoms.empty()) throw struct_parse_error("No atoms found in mmCIF _atom_site.");
}

// ------------------------------------------------------------------------------------------------
// Residue templates for heavy-atom-only residues. Only what changes typing is listed:
// D = heavy atom carrying polar hydrogen(s), A = nitrogen acceptor (oxygens are always OA),
// R = aromatic ring carbon (AutoDock type A).

enum : unsigned { TPL_DONOR = 1, TPL_ACCEPTOR = 2, TPL_AROMATIC = 4 };

struct residue_template {
    std::map<std::string, unsigned> atoms;
    bool amino_acid = false;
    bool ambiguous_his = false;  // HIS without a tautomer: ND1/NE2 may be donor or acceptor
    bool nucleotide = false;
};

typedef std::map<std::string, residue_template> template_map;

void add_templates(template_map& tm, const std::string& names, const std::string& spec,
                   bool amino_acid, bool backbone_donor = true, bool ambiguous_his = false) {
    residue_template tpl;
    tpl.amino_acid = amino_acid;
    tpl.ambiguous_his = ambiguous_his;
    if (amino_acid && backbone_donor) tpl.atoms["N"] = TPL_DONOR;

    std::istringstream atoms(spec);
    std::string item;
    while (atoms >> item) {
        const sz colon = item.find(':');
        unsigned flags = 0;
        for (char ch : item.substr(colon + 1)) {
            if (ch == 'D') flags |= TPL_DONOR;
            if (ch == 'A') flags |= TPL_ACCEPTOR;
            if (ch == 'R') flags |= TPL_AROMATIC;
        }
        tpl.atoms[item.substr(0, colon)] = flags;
    }
    std::istringstream res(names);
    std::string name;
    while (res >> name) tm[name] = tpl;
}

const template_map& residue_templates() {
    static const template_map tm = [] {
        template_map t;
        const std::string phe_ring = "CG:R CD1:R CD2:R CE1:R CE2:R CZ:R ";
        const std::string his_ring = "CG:R CD2:R CE1:R ";
        const std::string purine = "C2:R C4:R C5:R C6:R C8:R ";
        const std::string pyrimidine = "C2:R C4:R C5:R C6:R ";

        add_templates(t, "ALA GLY ILE LEU MET MSE VAL CYS CYX CYM ASP GLU", "", true);
        add_templates(t, "PRO", "", true, false);
        add_templates(t, "PHE", phe_ring, true);
        add_templates(t, "TYR", phe_ring + "OH:D", true);
        add_templates(t, "TRP", "CG:R CD1:R CD2:R CE2:R CE3:R CZ2:R CZ3:R CH2:R NE1:D", true);
        add_templates(t, "ARG", "NE:D NH1:D NH2:D", true);
        add_templates(t, "ASN", "ND2:D", true);
        add_templates(t, "GLN", "NE2:D", true);
        add_templates(t, "LYS LYN", "NZ:D", true);
        add_templates(t, "SER", "OG:D", true);
        add_templates(t, "THR", "OG1:D", true);
        add_templates(t, "ASH", "OD2:D", true);
        add_templates(t, "GLH", "OE2:D", true);
        add_templates(t, "HID HSD", his_ring + "ND1:D NE2:A", true);
        add_templates(t, "HIE HSE", his_ring + "ND1:A NE2:D", true);
        add_templates(t, "HIP HSP", his_ring + "ND1:D NE2:D", true);
        add_templates(t, "HIS", his_ring + "ND1:DA NE2:DA", true, true, true);

        add_templates(t, "HOH WAT DOD H2O TIP TIP3 SOL", "O:D", false);

        add_templates(t, "A", purine + "N1:A N3:A N7:A N6:D O2':D", false);
        add_templates(t, "G", purine + "N1:D N2:D N3:A N7:A O2':D", false);
        add_templates(t, "C", pyrimidine + "N3:A N4:D O2':D", false);
        add_templates(t, "U", pyrimidine + "N3:D O2':D", false);
        add_templates(t, "DA", purine + "N1:A N3:A N7:A N6:D", false);
        add_templates(t, "DG", purine + "N1:D N2:D N3:A N7:A", false);
        add_templates(t, "DC", pyrimidine + "N3:A N4:D", false);
        add_templates(t, "DT", pyrimidine + "N3:D", false);
        for (const char* nt : {"A", "C", "G", "U", "DA", "DC", "DG", "DT"}) t[nt].nucleotide = true;
        return t;
    }();
    return tm;
}

// ------------------------------------------------------------------------------------------------
// typing

// base AutoDock type from the element; false if AutoDock Vina cannot represent the element
bool element_to_atom_type(const std::string& el, atom& a) {
    static const char* const ad_elements[] = {"H",  "C",  "N",  "O",  "S",  "Se", "P",  "F",  "Cl",
                                              "Br", "I",  "Si", "At", "Mg", "Mn", "Zn", "Ca", "Fe"};
    a.ad = AD_TYPE_SIZE;
    if (el == "D") a.ad = AD_TYPE_H;
    for (const char* e : ad_elements)
        if (el == e) a.ad = string_to_ad_type(el);
    if (a.ad == AD_TYPE_O) a.ad = AD_TYPE_OA;
    if (is_non_ad_metal_name(el)) a.xs = XS_TYPE_Met_D;  // same as the PDBQT reader
    return a.acceptable_type();
}

bool is_metal(const atom& a) { return ad_type_to_el_type(a.ad) == EL_TYPE_Met || a.xs == XS_TYPE_Met_D; }

// uniform grid for neighbour searches; queries cover the 27 cells around a point
class spatial_grid {
public:
    spatial_grid(const atomv& atoms, fl cell) : m_cell(cell) {
        VINA_FOR_IN(i, atoms) m_cells[key_of(atoms[i].coords)].push_back(i);
    }

    template <typename F> void for_each_near(const vec& c, F f) const {
        const cell_key cx = index(c[0]), cy = index(c[1]), cz = index(c[2]);
        for (cell_key dx = -1; dx <= 1; ++dx)
            for (cell_key dy = -1; dy <= 1; ++dy)
                for (cell_key dz = -1; dz <= 1; ++dz) {
                    std::unordered_map<cell_key, szv>::const_iterator it
                        = m_cells.find(key(cx + dx, cy + dy, cz + dz));
                    if (it == m_cells.end()) continue;
                    for (sz j : it->second) f(j);
                }
    }

private:
    typedef long long cell_key;
    cell_key index(fl v) const { return cell_key(std::floor(v / m_cell)); }
    static cell_key key(cell_key x, cell_key y, cell_key z) {
        const cell_key off = 1 << 20;
        return ((x + off) << 42) | ((y + off) << 21) | (z + off);
    }
    cell_key key_of(const vec& c) const { return key(index(c[0]), index(c[1]), index(c[2])); }

    fl m_cell;
    std::unordered_map<cell_key, szv> m_cells;
};

struct connectivity {
    std::vector<int> h_parent;  // hydrogens: nearest bonded non-hydrogen, -1 if none
    szv heavy_degree;           // bonded non-hydrogen neighbours (metals included)
    szv h_count;                // hydrogens attached
    szv metal_neighbors;
    std::vector<szv> neighbors;  // bonded non-hydrogen, non-metal neighbours
};

// same bonding criterion as model::assign_bonds: r < 1.1 * (covalent radius i + covalent radius j)
connectivity perceive_bonds(const atomv& atoms, const spatial_grid& grid) {
    const sz n = atoms.size();
    connectivity con;
    con.h_parent.assign(n, -1);
    con.heavy_degree.assign(n, 0);
    con.h_count.assign(n, 0);
    con.metal_neighbors.assign(n, 0);
    con.neighbors.assign(n, szv());
    std::vector<fl> h_best(n, max_fl), h_best_metal(n, max_fl);
    std::vector<int> h_metal_parent(n, -1);

    VINA_FOR(i, n) {
        const atom& ai = atoms[i];
        grid.for_each_near(ai.coords, [&](sz j) {
            if (j <= i) return;
            const atom& aj = atoms[j];
            const bool hi = ai.is_hydrogen(), hj = aj.is_hydrogen();
            if (hi && hj) return;
            const fl r2 = vec_distance_sqr(ai.coords, aj.coords);
            if (r2 >= sqr(1.1 * ai.optimal_covalent_bond_length(aj))) return;
            if (hi || hj) {
                // a hydrogen pointing at a coordinated metal still belongs to its N/O, so
                // metals are only a fallback parent
                const sz h = hi ? i : j, heavy = hi ? j : i;
                const bool metal = is_metal(atoms[heavy]);
                std::vector<fl>& best = metal ? h_best_metal : h_best;
                std::vector<int>& parent = metal ? h_metal_parent : con.h_parent;
                if (r2 < best[h]) {
                    best[h] = r2;
                    parent[h] = int(heavy);
                }
            } else {
                ++con.heavy_degree[i];
                ++con.heavy_degree[j];
                const bool mi = is_metal(ai), mj = is_metal(aj);
                if (mj) ++con.metal_neighbors[i];
                if (mi) ++con.metal_neighbors[j];
                if (!mi && !mj) {
                    con.neighbors[i].push_back(j);
                    con.neighbors[j].push_back(i);
                }
            }
        });
    }
    VINA_FOR(i, n) {
        if (con.h_parent[i] < 0) con.h_parent[i] = h_metal_parent[i];
        if (con.h_parent[i] >= 0) ++con.h_count[sz(con.h_parent[i])];
    }
    return con;
}

// ------------------------------------------------------------------------------------------------
// Geometry-based perception of polar hydrogens for residues that have neither hydrogens nor a
// template (ligands, cofactors, modified residues). Bond lengths separate single from double
// bonds (C-OH >= 1.30 A > C=O, C#N <= 1.20 A), ring planarity identifies aromatic rings, and
// protonation follows pH 7 (carboxylates, phosphates and sulfates are deprotonated).

class geometry_typer {
public:
    geometry_typer(const atomv& atoms, const connectivity& con) : m_atoms(atoms), m_con(con) {}

    struct polar {
        bool donor = false;
        bool acceptor = false;
    };

    polar oxygen(sz i) const {
        polar t;
        t.acceptor = true;
        const szv& nb = m_con.neighbors[i];
        if (nb.empty()) {  // water or hydroxide
            t.donor = true;
            return t;
        }
        if (nb.size() > 1) return t;  // ether, ester, ring oxygen
        const sz x = nb[0];
        const sz xel = element(x);
        if (xel == EL_TYPE_C) {
            // carboxyl groups are carboxylates at pH 7; otherwise a long C-O bond is a hydroxyl
            t.donor = terminal_oxygens(x) == 1 && distance(i, x) > 1.30;
        } else if (xel == EL_TYPE_N) {
            t.donor = terminal_oxygens(x) == 1 && distance(i, x) >= 1.34;  // not nitro/N-oxide
        } else if (xel != EL_TYPE_P && xel != EL_TYPE_S) {
            t.donor = true;  // e.g. B-OH, Si-OH; phosphates and sulfates are deprotonated
        }
        return t;
    }

    polar nitrogen(sz i) const {
        polar t;
        if (m_con.metal_neighbors[i] > 0) return t;  // coordinated to a metal (e.g. heme)
        const szv& nb = m_con.neighbors[i];
        if (nb.empty()) {
            t.donor = true;  // ammonia/ammonium
        } else if (nb.size() == 1) {
            if (distance(i, nb[0]) <= 1.20) {
                t.acceptor = true;  // nitrile
            } else {
                t.donor = true;  // amine, amide or amidine NH2
                t.acceptor = is_sulfonyl(nb[0]);  // sulfonamide N, typed NA by Meeko
            }
        } else if (nb.size() == 2) {
            two_neighbor_nitrogen(i, t);
        } else if (nb.size() == 3) {
            t.donor = tertiary_aliphatic_amine(i);  // protonated at pH 7
        }
        // otherwise amide, aniline or aromatic N with three neighbours: plain N, as AutoDock
        return t;
    }

    bool aromatic_carbon(sz i) const {
        for (const szv& ring : planar_rings(i))
            if (ring.size() >= 5) return true;
        return false;
    }

private:
    sz element(sz i) const { return ad_type_to_el_type(m_atoms[i].ad); }
    fl distance(sz i, sz j) const {
        return std::sqrt(vec_distance_sqr(m_atoms[i].coords, m_atoms[j].coords));
    }

    // oxygens bonded only to atom x
    sz terminal_oxygens(sz x) const {
        sz count = 0;
        for (sz j : m_con.neighbors[x])
            if (element(j) == EL_TYPE_O && m_con.neighbors[j].size() == 1) ++count;
        return count;
    }

    bool is_carbonyl_carbon(sz c) const {
        if (element(c) != EL_TYPE_C) return false;
        for (sz j : m_con.neighbors[c])
            if (element(j) == EL_TYPE_O && m_con.neighbors[j].size() == 1 && distance(c, j) <= 1.30)
                return true;
        return false;
    }

    bool is_sulfonyl(sz x) const { return element(x) == EL_TYPE_S && terminal_oxygens(x) >= 2; }

    // carbon without double bonds, judged from the angles around it
    bool is_sp3_carbon(sz c) const {
        if (element(c) != EL_TYPE_C) return false;
        const szv& nb = m_con.neighbors[c];
        if (nb.size() == 2) return angle_deg(m_atoms[nb[0]].coords, m_atoms[c].coords, m_atoms[nb[1]].coords) < 117;
        if (nb.size() == 3) return angle_sum(c) < 350;
        return true;
    }

    fl angle_sum(sz center) const {
        const szv& nb = m_con.neighbors[center];
        fl sum = 0;
        VINA_FOR_IN(a, nb)
        VINA_RANGE(b, a + 1, nb.size())
        sum += angle_deg(m_atoms[nb[a]].coords, m_atoms[center].coords, m_atoms[nb[b]].coords);
        return sum;
    }

    // pyramidal N bonded by single bonds to three sp3 carbons
    bool tertiary_aliphatic_amine(sz n) const {
        if (angle_sum(n) >= 345) return false;
        for (sz c : m_con.neighbors[n])
            if (!is_sp3_carbon(c) || distance(n, c) < 1.43) return false;
        return true;
    }

    bool has_exocyclic_nitrogen(sz c, const szv& ring) const {
        for (sz j : m_con.neighbors[c])
            if (element(j) == EL_TYPE_N && std::find(ring.begin(), ring.end(), j) == ring.end())
                return true;
        return false;
    }

    // planar 5- or 6-membered rings through atom i
    std::vector<szv> planar_rings(sz start) const {
        std::vector<szv> rings;
        std::set<szv> seen;
        szv path(1, start);
        find_rings(start, path, rings, seen);
        return rings;
    }

    void find_rings(sz start, szv& path, std::vector<szv>& rings, std::set<szv>& seen) const {
        for (sz nb : m_con.neighbors[path.back()]) {
            if (nb == start && path.size() >= 5) {
                szv key = path;
                std::sort(key.begin(), key.end());
                if (seen.insert(key).second && planar(path)) rings.push_back(path);
                continue;
            }
            if (path.size() >= 6 || std::find(path.begin(), path.end(), nb) != path.end()) continue;
            if (m_con.neighbors[nb].size() > 3) continue;  // sp3 atom: not an aromatic ring
            path.push_back(nb);
            find_rings(start, path, rings, seen);
            path.pop_back();
        }
    }

    bool planar(const szv& ring) const {
        vec centroid(0, 0, 0);
        for (sz k : ring) centroid += m_atoms[k].coords;
        centroid = centroid / fl(ring.size());
        vec normal(0, 0, 0);
        VINA_FOR_IN(k, ring)
        normal += cross_product(m_atoms[ring[k]].coords - centroid,
                                m_atoms[ring[(k + 1) % ring.size()]].coords - centroid);
        const fl len = normal.norm();
        if (len < epsilon_fl) return false;
        normal = normal / len;
        for (sz k : ring)
            if (std::abs((m_atoms[k].coords - centroid) * normal) > 0.1) return false;
        return true;
    }

    void two_neighbor_nitrogen(sz i, polar& t) const {
        const szv& nb = m_con.neighbors[i];
        const std::vector<szv> rings = planar_rings(i);
        const szv* ring6 = nullptr;
        const szv* ring5 = nullptr;
        for (const szv& ring : rings) {
            if (ring.size() == 6 && !ring6) ring6 = &ring;
            if (ring.size() == 5 && !ring5) ring5 = &ring;
        }
        const sz carbonyls = sz(is_carbonyl_carbon(nb[0])) + sz(is_carbonyl_carbon(nb[1]));

        if (ring6) {
            // pyridine-like acceptor, unless next to a ring carbonyl (lactam N-H)
            const bool amino = has_exocyclic_nitrogen(nb[0], *ring6)
                               || has_exocyclic_nitrogen(nb[1], *ring6);
            t.acceptor = carbonyls == 0 || amino;
            t.donor = carbonyls > 0;
        } else if (ring5) {
            bool pyrrole_type_partner = false;
            sz pyridine_type_partners = 0;
            for (sz k : *ring5) {
                if (k == i) continue;
                const sz el = element(k);
                if (el == EL_TYPE_O || el == EL_TYPE_S
                    || (el == EL_TYPE_N && m_con.neighbors[k].size() == 3))
                    pyrrole_type_partner = true;
                else if (el == EL_TYPE_N && m_con.neighbors[k].size() == 2)
                    ++pyridine_type_partners;
            }
            if (carbonyls > 0) {
                t.donor = true;  // imide/lactam N-H
            } else if (pyrrole_type_partner) {
                t.acceptor = true;  // oxazole/thiazole N, N-substituted imidazole N3, purine N7
            } else if (pyridine_type_partners == 0) {
                t.donor = true;  // pyrrole/indole N-H
            } else {
                t.donor = t.acceptor = true;  // unsubstituted imidazole/pyrazole/triazole
            }
        } else if (angle_deg(m_atoms[nb[0]].coords, m_atoms[i].coords, m_atoms[nb[1]].coords)
                   >= 160) {
            t.acceptor = true;  // linear: azide, carbodiimide
        } else {
            // a short bond to carbon or nitrogen is a C=N / N=N double bond (no H) unless the
            // carbon is an amidine/guanidine, which is protonated at pH 7
            for (sz x : nb) {
                const sz xel = element(x);
                if ((xel == EL_TYPE_C && distance(i, x) <= 1.30)
                    || (xel == EL_TYPE_N && distance(i, x) <= 1.28)) {
                    bool amidine = false;
                    if (xel == EL_TYPE_C)
                        for (sz y : m_con.neighbors[x])
                            if (y != i && element(y) == EL_TYPE_N) amidine = true;
                    if (!amidine) {
                        t.acceptor = true;
                        return;
                    }
                }
            }
            t.donor = true;  // secondary amide, amine, sulfonamide, amidinium
            t.acceptor = is_sulfonyl(nb[0]) || is_sulfonyl(nb[1]);
        }
    }

    static fl angle_deg(const vec& a, const vec& center, const vec& b) {
        const vec u = a - center, v = b - center;
        fl c = (u * v) / (u.norm() * v.norm());
        c = std::max(fl(-1), std::min(fl(1), c));
        return std::acos(c) * 180 / pi;
    }

    const atomv& m_atoms;
    const connectivity& m_con;
};

// ------------------------------------------------------------------------------------------------
// Histidines without hydrogens: the tautomer is taken from hydrogen-bond partners in the
// direction where the ring N-H would point. A pure acceptor (carboxylate, backbone O) in front of
// a ring nitrogen means that nitrogen carries the hydrogen; a pure donor means it does not.
// Without evidence HID is used, the default of PDB2PQR and OpenMM/PDBFixer.

struct his_counts {
    sz hid = 0, hie = 0, hip = 0, unresolved = 0;
};

class histidine_resolver {
public:
    histidine_resolver(const atomv& atoms, const connectivity& con, const spatial_grid& grid,
                       const szv& residue_of, const std::vector<bool>& donor,
                       const std::vector<bool>& acceptor)
        : m_atoms(atoms),
          m_con(con),
          m_grid(grid),
          m_residue_of(residue_of),
          m_donor(donor),
          m_acceptor(acceptor) {}

    // +1: the nitrogen should carry the hydrogen, -1: it should not, 0: no evidence
    int vote(sz n) const {
        const szv& nb = m_con.neighbors[n];
        if (nb.size() != 2) return 0;
        const vec& p = m_atoms[n].coords;
        vec out = p - fl(0.5) * (m_atoms[nb[0]].coords + m_atoms[nb[1]].coords);
        const fl out_len = out.norm();
        if (out_len < epsilon_fl) return 0;
        out = out / out_len;

        fl best = 3.3;
        int result = 0;
        m_grid.for_each_near(p, [&](sz j) {
            if (m_residue_of[j] == m_residue_of[n]) return;
            const sz el = ad_type_to_el_type(m_atoms[j].ad);
            if (el != EL_TYPE_N && el != EL_TYPE_O) return;
            const fl d = std::sqrt(vec_distance_sqr(p, m_atoms[j].coords));
            if (d < 2.5 || d >= best) return;
            if (((m_atoms[j].coords - p) / d) * out < 0.64) return;  // more than ~50 deg off
            const bool jd = m_donor[j], ja = m_acceptor[j];
            if (jd == ja) return;  // water, hydroxyl, other histidine: no information
            best = d;
            result = ja ? 1 : -1;
        });
        return result;
    }

private:
    const atomv& m_atoms;
    const connectivity& m_con;
    const spatial_grid& m_grid;
    const szv& m_residue_of;
    const std::vector<bool>& m_donor;
    const std::vector<bool>& m_acceptor;
};

struct residue_info {
    std::string name;
    szv atoms;
    bool has_hydrogens = false;
};

std::string join_names(const std::set<std::string>& names) {
    std::string s;
    sz count = 0;
    for (const std::string& n : names) {
        if (count == 10) {
            s += ", ...";
            break;
        }
        if (count++) s += ", ";
        s += n;
    }
    return s;
}

}  // namespace

bool is_mmcif_file_name(const std::string& name) {
    return ends_with_ci(name, ".cif") || ends_with_ci(name, ".mmcif");
}

void parse_mmcif_rigid(const path& name, rigid& r) {
    ifile in(name);
    parse_mmcif_rigid(in, r);
}

void parse_mmcif_rigid(std::istream& in, rigid& r) {
    atom_site_reader reader;
    read_atom_site(in, reader);
    const std::vector<raw_atom>& raw = reader.atoms;
    const sz n = raw.size();

    atomv atoms(n);
    std::vector<residue_info> residues;
    szv residue_of(n);
    std::unordered_map<std::string, sz> residue_index;
    VINA_FOR(i, n) {
        atom& a = atoms[i];
        a.coords = raw[i].coords;
        a.charge = 0;
        a.number_sdf = 0;
        if (!element_to_atom_type(raw[i].element, a))
            throw struct_parse_error(
                "Element \"" + raw[i].element + "\" of atom " + raw[i].name + " in residue "
                + raw[i].res_name + " (mmCIF line " + std::to_string(raw[i].line)
                + ") cannot be represented by an AutoDock type.");

        std::unordered_map<std::string, sz>::iterator it = residue_index.find(raw[i].res_key);
        if (it == residue_index.end()) {
            it = residue_index.insert(std::make_pair(raw[i].res_key, residues.size())).first;
            residues.push_back(residue_info());
            residues.back().name = to_upper(raw[i].res_name);
        }
        residue_of[i] = it->second;
        residue_info& res = residues[it->second];
        res.atoms.push_back(i);
        if (a.is_hydrogen()) res.has_hydrogens = true;
    }

    const spatial_grid grid(atoms, 4.0);  // cell >= bond cutoff and H-bond search radius
    const connectivity con = perceive_bonds(atoms, grid);
    const geometry_typer geometry(atoms, con);
    const template_map& templates = residue_templates();
    std::vector<bool> donor(n, false), acceptor(n, false);
    std::set<std::string> geometry_residues;
    std::vector<const residue_info*> unresolved_his;

    for (const residue_info& res : residues) {
        template_map::const_iterator tpl_it = templates.find(res.name);
        const residue_template* tpl = tpl_it == templates.end() ? nullptr : &tpl_it->second;
        bool has_ca = false, has_c = false;
        for (sz i : res.atoms) {
            if (raw[i].name == "CA") has_ca = true;
            if (raw[i].name == "C") has_c = true;
        }
        const bool backbone_like = has_ca && has_c;

        for (sz i : res.atoms) {
            atom& a = atoms[i];
            const std::string& name = raw[i].name;
            const sz el = ad_type_to_el_type(a.ad);
            unsigned flags = 0;
            if (tpl) {
                std::map<std::string, unsigned>::const_iterator f = tpl->atoms.find(name);
                if (f != tpl->atoms.end()) flags = f->second;
            }

            if (el == EL_TYPE_H) {
                const int p = con.h_parent[i];
                const sz pel = p < 0 ? EL_TYPE_SIZE : ad_type_to_el_type(atoms[sz(p)].ad);
                a.ad = (pel == EL_TYPE_N || pel == EL_TYPE_O || pel == EL_TYPE_S) ? AD_TYPE_HD
                                                                                  : AD_TYPE_H;
            } else if (el == EL_TYPE_C) {
                if ((flags & TPL_AROMATIC) || (!tpl && geometry.aromatic_carbon(i)))
                    a.ad = AD_TYPE_A;
            } else if (el == EL_TYPE_S && a.ad == AD_TYPE_S && raw[i].element == "S") {
                if (con.heavy_degree[i] + con.h_count[i] < 4) a.ad = AD_TYPE_SA;
            } else if (el == EL_TYPE_N || el == EL_TYPE_O) {
                const bool backbone_n
                    = el == EL_TYPE_N && name == "N" && (backbone_like || (tpl && tpl->amino_acid));
                if (res.has_hydrogens) {
                    // prepared residue: donors come from the HD atoms, as for PDBQT. A
                    // nitrogen without H and with at most 2 covalent (non-metal) neighbours is
                    // an acceptor, except a backbone amide N (e.g. proline after a chain break)
                    donor[i] = con.h_count[i] > 0;
                    acceptor[i] = el == EL_TYPE_O
                                  || (con.h_count[i] == 0
                                      && con.heavy_degree[i] - con.metal_neighbors[i] <= 2
                                      && !backbone_n);
                } else if (tpl) {
                    donor[i] = (flags & TPL_DONOR) != 0;
                    acceptor[i] = el == EL_TYPE_O || (flags & TPL_ACCEPTOR);
                    if (el == EL_TYPE_N && con.metal_neighbors[i] > 0)
                        donor[i] = acceptor[i] = false;  // metal-coordinating nitrogen
                    // 5'/3' terminal hydroxyls of a nucleic acid chain
                    if (tpl->nucleotide && (name == "O5'" || name == "O3'")
                        && con.heavy_degree[i] <= 1)
                        donor[i] = true;
                } else {
                    const geometry_typer::polar p
                        = el == EL_TYPE_O ? geometry.oxygen(i) : geometry.nitrogen(i);
                    donor[i] = p.donor;
                    acceptor[i] = p.acceptor;
                    geometry_residues.insert(res.name);
                }
            }
        }

        if (tpl && tpl->ambiguous_his && !res.has_hydrogens) unresolved_his.push_back(&res);
    }

    // histidine tautomers, once every other residue has donors/acceptors
    his_counts his;
    const histidine_resolver resolver(atoms, con, grid, residue_of, donor, acceptor);
    std::vector<std::pair<sz, sz> > his_rings;
    std::vector<std::pair<int, int> > his_votes;
    for (const residue_info* res : unresolved_his) {
        int nd1 = -1, ne2 = -1;
        for (sz i : res->atoms) {
            if (raw[i].name == "ND1") nd1 = int(i);
            if (raw[i].name == "NE2") ne2 = int(i);
        }
        if (nd1 < 0 || ne2 < 0) continue;
        his_rings.push_back(std::make_pair(sz(nd1), sz(ne2)));
        his_votes.push_back(std::make_pair(resolver.vote(sz(nd1)), resolver.vote(sz(ne2))));
    }
    VINA_FOR_IN(k, his_rings) {
        const sz nd1 = his_rings[k].first, ne2 = his_rings[k].second;
        const bool nd1_metal = con.metal_neighbors[nd1] > 0, ne2_metal = con.metal_neighbors[ne2] > 0;
        int v1 = his_votes[k].first, v2 = his_votes[k].second;
        bool nd1_h, ne2_h;
        if (nd1_metal || ne2_metal) {
            // the ring nitrogen opposite a metal-bound one carries the hydrogen
            nd1_h = !nd1_metal;
            ne2_h = !ne2_metal;
        } else if (v1 > 0 && v2 > 0) {
            nd1_h = ne2_h = true;
            ++his.hip;
        } else if (v2 > 0 || v1 < 0) {
            nd1_h = false;
            ne2_h = true;
            ++his.hie;
        } else {
            // HID also when there is no evidence, as PDB2PQR and OpenMM/PDBFixer do
            nd1_h = true;
            ne2_h = false;
            if (v1 > 0 || v2 < 0)
                ++his.hid;
            else
                ++his.unresolved;
        }
        donor[nd1] = nd1_h;
        donor[ne2] = ne2_h;
        acceptor[nd1] = !nd1_h && !nd1_metal;
        acceptor[ne2] = !ne2_h && !ne2_metal;
    }

    VINA_FOR(i, n)
    if (acceptor[i] && ad_type_to_el_type(atoms[i].ad) == EL_TYPE_N) atoms[i].ad = AD_TYPE_NA;

    // donors are only needed where hydrogens are missing from the file
    VINA_FOR(i, n)
    if (residues[residue_of[i]].has_hydrogens) donor[i] = false;

    r.atoms.insert(r.atoms.end(), atoms.begin(), atoms.end());
    r.implicit_donors.resize(r.atoms.size() - n, false);
    r.implicit_donors.insert(r.implicit_donors.end(), donor.begin(), donor.end());

    if (!reader.skipped_models.empty())
        std::cerr << "WARNING: mmCIF receptor has " << reader.skipped_models.size() + 1
                  << " models; only the first one is used.\n";
    if (reader.skipped_altloc_atoms > 0)
        std::cerr << "WARNING: mmCIF receptor has alternate locations; kept the first one of each "
                     "residue (" << reader.skipped_altloc_atoms << " atoms ignored).\n";
    const sz n_his = his.hid + his.hie + his.hip + his.unresolved;
    if (n_his > 0)
        std::cerr << "NOTE: " << n_his
                  << " histidines without hydrogens; tautomers assigned from hydrogen-bond "
                     "partners (HID " << his.hid << ", HIE " << his.hie << ", HIP " << his.hip
                  << ", HID by default " << his.unresolved << ").\n";
    if (!geometry_residues.empty())
        std::cerr << "WARNING: mmCIF receptor residues without hydrogens and without a built-in "
                     "template (" << join_names(geometry_residues)
                  << ") were typed from their geometry (bond lengths, angles, ring planarity). "
                     "Add hydrogens to the receptor for exact typing.\n";
}
