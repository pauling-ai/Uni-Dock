#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <fstream>
#include <sstream>
#include "parse_pdbqt.h"
#include "parse_mmcif.h"
#include "parse_error.h"
#include "atom.h"

namespace {

struct pdbqt_atom_label {
    std::string name, res_name;
    vec coords;
};

// residue/atom names of def.pdbqt, to explain typing differences by atom
std::vector<pdbqt_atom_label> read_pdbqt_labels(const std::string& file) {
    std::vector<pdbqt_atom_label> labels;
    std::ifstream in(file);
    std::string line;
    while (std::getline(in, line)) {
        if (line.compare(0, 4, "ATOM") != 0 && line.compare(0, 6, "HETATM") != 0) continue;
        pdbqt_atom_label l;
        std::istringstream(line.substr(12, 4)) >> l.name;
        std::istringstream(line.substr(17, 3)) >> l.res_name;
        l.coords = vec(std::stof(line.substr(30, 8)), std::stof(line.substr(38, 8)),
                       std::stof(line.substr(46, 8)));
        labels.push_back(l);
    }
    return labels;
}

bool same_coords(const vec& a, const vec& b) { return vec_distance_sqr(a, b) < 1e-6; }

}  // namespace

TEST_CASE("mmcif file name detection", "[parse_mmcif]") {
    REQUIRE(is_mmcif_file_name("receptor.cif"));
    REQUIRE(is_mmcif_file_name("dir.v2/receptor.CIF"));
    REQUIRE(is_mmcif_file_name("receptor.mmcif"));
    REQUIRE_FALSE(is_mmcif_file_name("receptor.pdbqt"));
    REQUIRE_FALSE(is_mmcif_file_name("receptor.pdb"));
    REQUIRE_FALSE(is_mmcif_file_name(""));
}

TEST_CASE("parse mmcif rigid", "[parse_mmcif_rigid]") {
    rigid r;
    REQUIRE_NOTHROW(parse_mmcif_rigid(path("./test_data/def.cif"), r));
    REQUIRE(r.atoms.size() == 1613);
    REQUIRE(r.implicit_donors.size() == r.atoms.size());
}

TEST_CASE("mmcif receptor with hydrogens is typed exactly like the pdbqt receptor",
          "[parse_receptor_mmcif]") {
    model from_pdbqt = parse_receptor_pdbqt("./test_data/def.pdbqt");
    model from_cif = parse_receptor_mmcif("./test_data/def.cif");

    REQUIRE(from_cif.grid_atoms.size() == from_pdbqt.grid_atoms.size());
    VINA_FOR_IN(i, from_cif.grid_atoms) {
        const atom& c = from_cif.grid_atoms[i];
        const atom& p = from_pdbqt.grid_atoms[i];
        INFO("atom index " << i);
        REQUIRE(same_coords(c.coords, p.coords));
        REQUIRE(c.el == p.el);
        REQUIRE(c.xs == p.xs);
        REQUIRE(c.bonds.size() == p.bonds.size());
    }
}

TEST_CASE("mmcif receptor without hydrogens uses residue templates",
          "[parse_receptor_mmcif]") {
    model from_pdbqt = parse_receptor_pdbqt("./test_data/def.pdbqt");
    model from_cif = parse_receptor_mmcif("./test_data/def_noH.cif");
    const std::vector<pdbqt_atom_label> labels = read_pdbqt_labels("./test_data/def.pdbqt");
    REQUIRE(labels.size() == from_pdbqt.grid_atoms.size());

    sz heavy = 0, cif_donors = 0, pdbqt_donors = 0;
    std::vector<std::string> mismatches;
    VINA_FOR_IN(i, from_pdbqt.grid_atoms) {
        const atom& p = from_pdbqt.grid_atoms[i];
        if (p.is_hydrogen()) continue;
        REQUIRE(heavy < from_cif.grid_atoms.size());
        const atom& c = from_cif.grid_atoms[heavy++];
        REQUIRE(same_coords(c.coords, p.coords));
        if (xs_is_donor(c.xs)) ++cif_donors;
        if (xs_is_donor(p.xs)) ++pdbqt_donors;
        if (c.xs != p.xs) mismatches.push_back(labels[i].res_name + ":" + labels[i].name);
    }
    REQUIRE(heavy == from_cif.grid_atoms.size());

    // The only atoms a template cannot type like the prepared PDBQT are histidine ring
    // nitrogens: HIS has no tautomer information (donor+acceptor instead of donor), and HIZ is
    // not a standard residue name (its 2 residues lose the ring donors).
    for (const std::string& m : mismatches) {
        INFO(m);
        const bool his_ring_n = m == "HIS:ND1" || m == "HIS:NE2" || m == "HIZ:ND1" || m == "HIZ:NE2";
        CHECK(his_ring_n);
    }
    CHECK(mismatches.size() == 8);
    CHECK(cif_donors == pdbqt_donors - 4);
}

TEST_CASE("mmcif tokenizer, models and alternate locations", "[parse_mmcif_rigid]") {
    const std::string cif = "data_t\n"
                            "# comment\n"
                            "_struct.title\n"
                            ";multi-line\n"
                            "text with 'quotes' and _tags\n"
                            ";\n"
                            "loop_\n"
                            "_atom_site.type_symbol\n"
                            "_atom_site.label_atom_id\n"
                            "_atom_site.label_alt_id\n"
                            "_atom_site.label_comp_id\n"
                            "_atom_site.auth_seq_id\n"
                            "_atom_site.Cartn_x\n"
                            "_atom_site.Cartn_y\n"
                            "_atom_site.Cartn_z\n"
                            "_atom_site.pdbx_PDB_model_num\n"
                            "O O     . HOH 1 0.0  0 0 1\n"
                            "O O     A HOH 2 5.0  0 0 1\n"
                            "O O     B HOH 2 5.5  0 0 1\n"
                            "N 'N'   . GLY 3 10.0 0 0 1\n"
                            "C \"C1'\" . DA  4 20.0 0 0 1\n"
                            "O O     . HOH 1 0.0  0 0 2\n"
                            "#\n"
                            "loop_\n"
                            "_other.value\n"
                            "1\n";
    std::istringstream in(cif);
    rigid r;
    REQUIRE_NOTHROW(parse_mmcif_rigid(in, r));
    REQUIRE(r.atoms.size() == 4);
    REQUIRE(r.atoms[1].coords[0] == Catch::Approx(5.0).margin(1e-6));
    REQUIRE(r.atoms[0].ad == AD_TYPE_OA);
    REQUIRE(r.atoms[2].ad == AD_TYPE_N);
    REQUIRE(r.atoms[3].ad == AD_TYPE_C);
    REQUIRE(r.implicit_donors == std::vector<bool>({true, true, true, false}));
}

TEST_CASE("mmcif single atom_site without loop and metals", "[parse_mmcif_rigid]") {
    std::istringstream in("data_zn\n"
                          "_atom_site.type_symbol ZN\n"
                          "_atom_site.label_atom_id ZN\n"
                          "_atom_site.label_comp_id ZN\n"
                          "_atom_site.Cartn_x 1.0\n"
                          "_atom_site.Cartn_y 2.0\n"
                          "_atom_site.Cartn_z 3.0(2)\n");
    rigid r;
    REQUIRE_NOTHROW(parse_mmcif_rigid(in, r));
    REQUIRE(r.atoms.size() == 1);
    REQUIRE(r.atoms[0].ad == AD_TYPE_Zn);
    REQUIRE(r.atoms[0].coords[2] == Catch::Approx(3.0).margin(1e-6));
}

TEST_CASE("mmcif errors", "[parse_mmcif_rigid]") {
    rigid r;
    std::istringstream no_atoms("data_x\n_entry.id X\n");
    REQUIRE_THROWS_AS(parse_mmcif_rigid(no_atoms, r), struct_parse_error);

    std::istringstream bad_element("data_x\nloop_\n_atom_site.type_symbol\n"
                                   "_atom_site.label_atom_id\n_atom_site.Cartn_x\n"
                                   "_atom_site.Cartn_y\n_atom_site.Cartn_z\nX X1 0 0 0\n");
    REQUIRE_THROWS_AS(parse_mmcif_rigid(bad_element, r), struct_parse_error);

    std::istringstream ragged("data_x\nloop_\n_atom_site.type_symbol\n_atom_site.label_atom_id\n"
                              "_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\n"
                              "C C1 0 0\n");
    REQUIRE_THROWS_AS(parse_mmcif_rigid(ragged, r), struct_parse_error);
}
