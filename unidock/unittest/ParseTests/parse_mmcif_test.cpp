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

// "type_symbol atom_id comp_id seq x y z" rows -> mmCIF text
std::string atom_site_cif(const std::vector<std::string>& rows) {
    std::string cif = "data_t\nloop_\n_atom_site.type_symbol\n_atom_site.label_atom_id\n"
                      "_atom_site.label_comp_id\n_atom_site.auth_seq_id\n_atom_site.Cartn_x\n"
                      "_atom_site.Cartn_y\n_atom_site.Cartn_z\n";
    for (const std::string& row : rows) cif += row + "\n";
    return cif;
}

rigid parse_rows(const std::vector<std::string>& rows) {
    std::istringstream in(atom_site_cif(rows));
    rigid r;
    parse_mmcif_rigid(in, r);
    REQUIRE(r.atoms.size() == rows.size());
    REQUIRE(r.implicit_donors.size() == rows.size());
    return r;
}

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

// Atoms in the template tests are 10 A apart, so no bonds are perceived between them.
TEST_CASE("mmcif amino acid templates without hydrogens", "[parse_mmcif_templates]") {
    rigid r = parse_rows({
        "N N   ALA 1 0 0 0",     //  0 backbone N: donor
        "O O   ALA 1 10 0 0",    //  1
        "N N   PRO 2 20 0 0",    //  2 proline N: no H
        "N ND1 HID 3 30 0 0",    //  3
        "N NE2 HID 3 40 0 0",    //  4
        "N ND1 HIE 4 50 0 0",    //  5
        "N NE2 HIE 4 60 0 0",    //  6
        "N ND1 HIP 5 70 0 0",    //  7
        "N NE2 HIP 5 80 0 0",    //  8
        "N ND1 HIS 6 90 0 0",    //  9
        "N NE2 HIS 6 100 0 0",   // 10
        "O OD2 ASH 7 110 0 0",   // 11
        "O OD2 ASP 8 120 0 0",   // 12
        "O OE2 GLH 9 130 0 0",   // 13
        "O OG  SER 10 140 0 0",  // 14
        "N NZ  LYS 11 150 0 0",  // 15
        "C CZ  PHE 12 160 0 0",  // 16 aromatic
        "C CB  PHE 12 170 0 0",  // 17
        "N NE1 TRP 13 180 0 0",  // 18
        "SE SE MSE 14 190 0 0",  // 19 selenium behaves as sulfur
        "O O   HOH 15 200 0 0",   // 20 water
        "N N   UNK 16 210 0 0",   // 21 unknown residue, no backbone
    });
    const std::vector<bool> expected_donors = {true,  false, false, true,  false, false,
                                               true,  true,  true,  true,  true,  true,
                                               false, true,  true,  true,  false, false,
                                               true,  false, true,  false};
    REQUIRE(r.implicit_donors == expected_donors);

    REQUIRE(r.atoms[0].ad == AD_TYPE_N);
    REQUIRE(r.atoms[1].ad == AD_TYPE_OA);
    REQUIRE(r.atoms[3].ad == AD_TYPE_N);   // HID ND1 protonated
    REQUIRE(r.atoms[4].ad == AD_TYPE_NA);  // HID NE2 acceptor
    REQUIRE(r.atoms[5].ad == AD_TYPE_NA);  // HIE ND1 acceptor
    REQUIRE(r.atoms[6].ad == AD_TYPE_N);
    REQUIRE(r.atoms[7].ad == AD_TYPE_N);
    REQUIRE(r.atoms[8].ad == AD_TYPE_N);
    REQUIRE(r.atoms[9].ad == AD_TYPE_NA);  // HIS: donor and acceptor
    REQUIRE(r.atoms[10].ad == AD_TYPE_NA);
    REQUIRE(r.atoms[16].ad == AD_TYPE_A);
    REQUIRE(r.atoms[17].ad == AD_TYPE_C);
    REQUIRE(r.atoms[19].ad == AD_TYPE_S);
    REQUIRE(r.atoms[21].ad == AD_TYPE_N);
}

TEST_CASE("mmcif nucleotide templates without hydrogens", "[parse_mmcif_templates]") {
    rigid r = parse_rows({
        "N N1    A  1 0 0 0",    // 0 acceptor
        "N N6    A  1 10 0 0",   // 1 donor
        "C C8    A  1 20 0 0",   // 2 aromatic
        "O \"O2'\" A  1 30 0 0", // 3 RNA 2'-OH donor
        "N N1    G  2 40 0 0",   // 4 donor
        "N N7    G  2 50 0 0",   // 5 acceptor
        "N N3    C  3 60 0 0",   // 6 acceptor
        "N N4    C  3 70 0 0",   // 7 donor
        "N N3    U  4 80 0 0",   // 8 donor
        "N N3    DT 5 90 0 0",   // 9 donor
        "N N3    DA 6 100 0 0",  // 10 acceptor
        "O O4*   DT 5 110 0 0",  // 11 old-style primed name
        "O OP1   DA 6 120 0 0",  // 12
    });
    const std::vector<bool> expected_donors = {false, true,  false, true,  true,  false, false,
                                               true,  true,  true,  false, false, false};
    REQUIRE(r.implicit_donors == expected_donors);
    REQUIRE(r.atoms[0].ad == AD_TYPE_NA);
    REQUIRE(r.atoms[1].ad == AD_TYPE_N);
    REQUIRE(r.atoms[2].ad == AD_TYPE_A);
    REQUIRE(r.atoms[5].ad == AD_TYPE_NA);
    REQUIRE(r.atoms[6].ad == AD_TYPE_NA);
    REQUIRE(r.atoms[10].ad == AD_TYPE_NA);
    REQUIRE(r.atoms[11].ad == AD_TYPE_OA);
}

TEST_CASE("mmcif histidine bound to a metal", "[parse_mmcif_templates]") {
    // ZN 2.1 A from NE2: NE2 is neither donor nor acceptor, ND1 carries the hydrogen
    rigid his = parse_rows({
        "N ND1 HIS 1 0 0 0",
        "N NE2 HIS 1 10 0 0",
        "ZN ZN ZN 2 12.1 0 0",
    });
    REQUIRE(his.implicit_donors == std::vector<bool>({true, false, false}));
    REQUIRE(his.atoms[0].ad == AD_TYPE_N);
    REQUIRE(his.atoms[1].ad == AD_TYPE_N);
    REQUIRE(his.atoms[2].ad == AD_TYPE_Zn);

    // a named tautomer keeps its template, but the metal-bound nitrogen loses donor/acceptor
    rigid hie = parse_rows({
        "N ND1 HIE 1 0 0 0",
        "N NE2 HIE 1 10 0 0",
        "ZN ZN ZN 2 -2.1 0 0",
    });
    REQUIRE(hie.implicit_donors == std::vector<bool>({false, true, false}));
    REQUIRE(hie.atoms[0].ad == AD_TYPE_N);

    // non-AutoDock metals are accepted like in PDBQT
    rigid cu = parse_rows({"CU CU CU 1 0 0 0"});
    REQUIRE(cu.atoms[0].xs == XS_TYPE_Met_D);
}

TEST_CASE("mmcif residues with hydrogens are typed from connectivity", "[parse_mmcif_rigid]") {
    rigid r = parse_rows({
        "C C1  LIG 1 0.000 0.000 0.000",
        "N N1  LIG 1 1.330 0.000 0.000",   // 2 heavy neighbours, no H -> NA
        "C C2  LIG 1 2.000 1.150 0.000",
        "N N2  LIG 1 -0.700 -1.150 0.000", // has H -> N (donor through HD)
        "H H2  LIG 1 -0.200 -2.000 0.000", // on N -> HD
        "H HC  LIG 1 -0.500 0.900 0.000",  // on C -> H
        "O O1  LIG 1 3.300 1.150 0.000",
    });
    REQUIRE(r.atoms[0].ad == AD_TYPE_C);
    REQUIRE(r.atoms[1].ad == AD_TYPE_NA);
    REQUIRE(r.atoms[3].ad == AD_TYPE_N);
    REQUIRE(r.atoms[4].ad == AD_TYPE_HD);
    REQUIRE(r.atoms[5].ad == AD_TYPE_H);
    REQUIRE(r.atoms[6].ad == AD_TYPE_OA);
    // donors come from HD atoms, never from templates, when hydrogens are present
    REQUIRE(r.implicit_donors == std::vector<bool>(7, false));
}

TEST_CASE("mmcif without type_symbol or label columns", "[parse_mmcif_rigid]") {
    std::istringstream in("data_t\nloop_\n_atom_site.auth_atom_id\n_atom_site.auth_comp_id\n"
                          "_atom_site.auth_asym_id\n_atom_site.auth_seq_id\n"
                          "_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\n"
                          "CA ALA A 1 0 0 0\n"
                          "CA CA  A 2 10 0 0\n"
                          "ZN ZN  A 3 20 0 0\n"
                          "1HB ALA A 1 30 0 0\n");
    rigid r;
    REQUIRE_NOTHROW(parse_mmcif_rigid(in, r));
    REQUIRE(r.atoms.size() == 4);
    REQUIRE(r.atoms[0].ad == AD_TYPE_C);   // alpha carbon
    REQUIRE(r.atoms[1].ad == AD_TYPE_Ca);  // calcium ion
    REQUIRE(r.atoms[2].ad == AD_TYPE_Zn);
    REQUIRE(r.atoms[3].ad == AD_TYPE_H);
}
