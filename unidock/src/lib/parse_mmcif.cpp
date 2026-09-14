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
     H bonded to N/O/S -> HD, N without H and with <= 2 neighbours -> NA, O -> OA, S -> SA.
     Donors then come out of model::assign_types through the HD atoms, exactly as for PDBQT.
   - residues without hydrogens (typical for PDB/AlphaFold files) use built-in templates of
     standard amino acids, nucleotides and water to know which heavy atoms carry polar
     hydrogens. Those atoms are reported in rigid::implicit_donors and promoted to XS donors
     after typing (see parse_receptor_mmcif).

*/

#include "parse_mmcif.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdlib>
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

struct connectivity {
    std::vector<int> h_parent;  // hydrogens: nearest bonded non-hydrogen, -1 if none
    szv heavy_degree;           // bonded non-hydrogen neighbours (metals included)
    szv h_count;                // hydrogens attached
    szv metal_neighbors;
};

// same bonding criterion as model::assign_bonds: r < 1.1 * (covalent radius i + covalent radius j)
connectivity perceive_bonds(const atomv& atoms) {
    const sz n = atoms.size();
    connectivity con;
    con.h_parent.assign(n, -1);
    con.heavy_degree.assign(n, 0);
    con.h_count.assign(n, 0);
    con.metal_neighbors.assign(n, 0);
    std::vector<fl> h_best(n, max_fl);

    const fl cell = 4.0;  // >= 1.1 * 2 * largest covalent radius
    typedef long long cell_key;
    const auto cell_index = [&](fl c) { return cell_key(std::floor(c / cell)); };
    const auto key = [](cell_key x, cell_key y, cell_key z) {
        const cell_key off = 1 << 20;
        return ((x + off) << 42) | ((y + off) << 21) | (z + off);
    };
    std::unordered_map<cell_key, szv> cells;
    VINA_FOR(i, n) {
        const vec& c = atoms[i].coords;
        cells[key(cell_index(c[0]), cell_index(c[1]), cell_index(c[2]))].push_back(i);
    }

    VINA_FOR(i, n) {
        const atom& ai = atoms[i];
        const vec& c = ai.coords;
        const cell_key cx = cell_index(c[0]), cy = cell_index(c[1]), cz = cell_index(c[2]);
        for (cell_key dx = -1; dx <= 1; ++dx)
            for (cell_key dy = -1; dy <= 1; ++dy)
                for (cell_key dz = -1; dz <= 1; ++dz) {
                    std::unordered_map<cell_key, szv>::const_iterator it
                        = cells.find(key(cx + dx, cy + dy, cz + dz));
                    if (it == cells.end()) continue;
                    for (sz j : it->second) {
                        if (j <= i) continue;
                        const atom& aj = atoms[j];
                        const bool hi = ai.is_hydrogen(), hj = aj.is_hydrogen();
                        if (hi && hj) continue;
                        const fl r2 = vec_distance_sqr(c, aj.coords);
                        if (r2 >= sqr(1.1 * ai.optimal_covalent_bond_length(aj))) continue;
                        if (hi || hj) {
                            const sz h = hi ? i : j, heavy = hi ? j : i;
                            if (r2 < h_best[h]) {
                                h_best[h] = r2;
                                con.h_parent[h] = int(heavy);
                            }
                        } else {
                            ++con.heavy_degree[i];
                            ++con.heavy_degree[j];
                            if (is_metal(aj)) ++con.metal_neighbors[i];
                            if (is_metal(ai)) ++con.metal_neighbors[j];
                        }
                    }
                }
    }
    VINA_FOR(i, n)
    if (con.h_parent[i] >= 0) ++con.h_count[sz(con.h_parent[i])];
    return con;
}

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
        residue_info& res = residues[it->second];
        res.atoms.push_back(i);
        if (a.is_hydrogen()) res.has_hydrogens = true;
    }

    const connectivity con = perceive_bonds(atoms);
    const template_map& templates = residue_templates();
    std::vector<bool> donor(n, false), acceptor(n, false);
    std::set<std::string> untyped_residues, ambiguous_his;

    for (const residue_info& res : residues) {
        template_map::const_iterator tpl_it = templates.find(res.name);
        const residue_template* tpl = tpl_it == templates.end() ? nullptr : &tpl_it->second;
        bool backbone_like = false, has_ca = false, has_c = false;
        for (sz i : res.atoms) {
            if (raw[i].name == "CA") has_ca = true;
            if (raw[i].name == "C") has_c = true;
        }
        backbone_like = has_ca && has_c;

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
                if (flags & TPL_AROMATIC) a.ad = AD_TYPE_A;
            } else if (el == EL_TYPE_S && a.ad == AD_TYPE_S && raw[i].element == "S") {
                if (con.heavy_degree[i] + con.h_count[i] < 4) a.ad = AD_TYPE_SA;
            } else if (el == EL_TYPE_N || el == EL_TYPE_O) {
                if (res.has_hydrogens) {
                    // prepared residue: donors come from the HD atoms, as for PDBQT
                    if (el == EL_TYPE_N)
                        acceptor[i] = con.h_count[i] == 0 && con.heavy_degree[i] <= 2;
                } else if (tpl) {
                    donor[i] = (flags & TPL_DONOR) != 0;
                    acceptor[i] = el == EL_TYPE_N && (flags & TPL_ACCEPTOR);
                    if (el == EL_TYPE_N && con.metal_neighbors[i] > 0)
                        donor[i] = acceptor[i] = false;  // metal-coordinating nitrogen
                } else {
                    if (el == EL_TYPE_N && name == "N" && backbone_like
                        && con.heavy_degree[i] <= 2)
                        donor[i] = true;  // backbone of a modified amino acid
                    else if (!(el == EL_TYPE_O && (name == "O" || name == "OXT") && backbone_like))
                        untyped_residues.insert(res.name);
                }
            }
        }

        if (tpl && tpl->ambiguous_his && !res.has_hydrogens) {
            int nd1 = -1, ne2 = -1;
            for (sz i : res.atoms) {
                if (raw[i].name == "ND1") nd1 = int(i);
                if (raw[i].name == "NE2") ne2 = int(i);
            }
            // the ring nitrogen opposite a metal-bound one carries the hydrogen
            if (nd1 >= 0 && ne2 >= 0 && con.metal_neighbors[sz(nd1)] > 0) {
                donor[sz(ne2)] = true;
                acceptor[sz(ne2)] = false;
            } else if (nd1 >= 0 && ne2 >= 0 && con.metal_neighbors[sz(ne2)] > 0) {
                donor[sz(nd1)] = true;
                acceptor[sz(nd1)] = false;
            } else {
                ambiguous_his.insert(res.name);
            }
        }
    }

    VINA_FOR(i, n)
    if (acceptor[i]) atoms[i].ad = AD_TYPE_NA;

    r.atoms.insert(r.atoms.end(), atoms.begin(), atoms.end());
    r.implicit_donors.resize(r.atoms.size() - n, false);
    r.implicit_donors.insert(r.implicit_donors.end(), donor.begin(), donor.end());

    if (!reader.skipped_models.empty())
        std::cerr << "WARNING: mmCIF receptor has " << reader.skipped_models.size() + 1
                  << " models; only the first one is used.\n";
    if (reader.skipped_altloc_atoms > 0)
        std::cerr << "WARNING: mmCIF receptor has alternate locations; kept the first one of each "
                     "residue (" << reader.skipped_altloc_atoms << " atoms ignored).\n";
    if (!ambiguous_his.empty())
        std::cerr << "WARNING: mmCIF receptor has histidines without hydrogens; ND1 and NE2 are "
                     "treated as both donor and acceptor. Add hydrogens or name them HID/HIE/HIP "
                     "to set the tautomer.\n";
    if (!untyped_residues.empty())
        std::cerr << "WARNING: mmCIF receptor residues without hydrogens and without a built-in "
                     "template (" << join_names(untyped_residues)
                  << "): their hydrogen-bond donors cannot be identified. Add hydrogens to the "
                     "receptor for exact typing.\n";
}
