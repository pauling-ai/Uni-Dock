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

    // The only atoms typed differently from this PDBQT are histidine ring nitrogens: the PDBQT
    // protonates both ring nitrogens of every histidine, while without hydrogens each HIS gets
    // one tautomer from its hydrogen-bond partners, and the two HIZ
    // (not a standard name) are typed from geometry as imidazoles (donor and acceptor).
    for (const std::string& m : mismatches) {
        INFO(m);
        const bool his_ring_n = m == "HIS:ND1" || m == "HIS:NE2" || m == "HIZ:ND1" || m == "HIZ:NE2";
        CHECK(his_ring_n);
    }
    CHECK(mismatches.size() == 6);
    CHECK(cif_donors == pdbqt_donors - 2);
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
        "N N   UNK 16 210 0 0",   // 21 no template, no neighbours: ammonia-like donor
    });
    const std::vector<bool> expected_donors = {true,  false, false, true,  false, false,
                                               true,  true,  true,  true,  false, true,
                                               false, true,  true,  true,  false, false,
                                               true,  false, true,  true};
    REQUIRE(r.implicit_donors == expected_donors);

    REQUIRE(r.atoms[0].ad == AD_TYPE_N);
    REQUIRE(r.atoms[1].ad == AD_TYPE_OA);
    REQUIRE(r.atoms[3].ad == AD_TYPE_N);   // HID ND1 protonated
    REQUIRE(r.atoms[4].ad == AD_TYPE_NA);  // HID NE2 acceptor
    REQUIRE(r.atoms[5].ad == AD_TYPE_NA);  // HIE ND1 acceptor
    REQUIRE(r.atoms[6].ad == AD_TYPE_N);
    REQUIRE(r.atoms[7].ad == AD_TYPE_N);
    REQUIRE(r.atoms[8].ad == AD_TYPE_N);
    REQUIRE(r.atoms[9].ad == AD_TYPE_N);  // HIS without H-bond partners: HID
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

TEST_CASE("mmcif hydrogen typing next to metals and chain breaks", "[parse_mmcif_rigid]") {
    SECTION("backbone N without H is not an acceptor (proline after a chain break)") {
        rigid r = parse_rows({
            "N N  PRO 1 0.000 0.000 0.000",
            "C CA PRO 1 1.470 0.000 0.000",
            "C CD PRO 1 -0.900 1.100 0.000",
            "H HA PRO 1 1.800 1.000 0.000",
        });
        REQUIRE(r.atoms[0].ad == AD_TYPE_N);
        REQUIRE(r.atoms[3].ad == AD_TYPE_H);
    }
    SECTION("hydrogen pointing at a coordinated metal stays polar") {
        rigid r = parse_rows({
            "N  NE2 HIE 1 0.000 0.000 0.000",
            "C  CE1 HIE 1 -0.700 1.100 0.000",
            "C  CD2 HIE 1 -0.700 -1.100 0.000",
            "H  HE2 HIE 1 1.020 0.000 0.000",   // 0.97 A from ZN, 1.02 A from NE2
            "ZN ZN  ZN  2 1.990 0.000 0.000",
        });
        REQUIRE(r.atoms[3].ad == AD_TYPE_HD);
        REQUIRE(r.atoms[0].ad == AD_TYPE_N);
    }
    SECTION("metal-coordinated N without H is still an acceptor, like Meeko") {
        rigid r = parse_rows({
            "N  N1 LIG 1 0.000 0.000 0.000",
            "C  C1 LIG 1 1.300 0.000 0.000",
            "C  C2 LIG 1 -0.650 1.130 0.000",
            "H  H1 LIG 1 1.800 0.900 0.000",
            "ZN ZN ZN  2 -0.650 -1.900 0.000",
        });
        REQUIRE(r.atoms[0].ad == AD_TYPE_NA);
        REQUIRE(r.atoms[3].ad == AD_TYPE_H);
    }
    SECTION("nucleic acid terminal hydroxyls without hydrogens are donors") {
        rigid r = parse_rows({
            "O \"O5'\" DC 1 0.000 0.000 0.000",   // 5' end: only C5'
            "C \"C5'\" DC 1 1.430 0.000 0.000",
            "O \"O5'\" DC 2 20.000 0.000 0.000",  // linked to the phosphate
            "C \"C5'\" DC 2 21.430 0.000 0.000",
            "P P       DC 2 18.400 0.000 0.000",
        });
        REQUIRE(r.implicit_donors == std::vector<bool>({true, false, false, false, false}));
    }
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

namespace {

// index of the atom "res seq name" in rows built by atom_site_cif
sz row_index(const std::vector<std::string>& rows, const std::string& res, int seq,
             const std::string& name) {
    VINA_FOR_IN(i, rows) {
        std::istringstream fields(rows[i]);
        std::string el, atom_name, res_name;
        int res_seq;
        fields >> el >> atom_name >> res_name >> res_seq;
        if (atom_name == name && res_name == res && res_seq == seq) return i;
    }
    FAIL("atom " << res << " " << seq << " " << name << " not found");
    return 0;
}

}  // namespace

TEST_CASE("mmcif residues without template or hydrogens are typed from geometry",
          "[parse_mmcif_geometry]") {
    // heavy atoms of RDKit/MMFF-optimised molecules, one residue each
    const std::vector<std::string> rows = {
        "C C1 ETH 1 -0.863 0.220 0.203",
        "C C2 ETH 1 0.406 -0.375 -0.370",
        "O O1 ETH 1 1.527 0.404 0.022",
        "C C1 ACN 2 23.733 -0.034 -0.251",
        "C C2 ACN 2 24.960 -0.009 0.620",
        "C C3 ACN 2 26.289 0.039 -0.084",
        "O O1 ACN 2 24.880 -0.028 1.847",
        "C C1 ACT 3 49.362 -0.055 -0.005",
        "C C2 ACT 3 50.875 0.069 0.007",
        "O O1 ACT 3 51.507 -1.018 -0.108",
        "O O2 ACT 3 51.309 1.249 0.130",
        "C C1 ACI 4 74.043 -0.073 -0.088",
        "C C2 ACI 4 75.490 0.288 -0.018",
        "O O1 ACI 4 75.976 1.384 -0.239",
        "O O2 ACI 4 76.264 -0.756 0.331",
        "C C1 PHN 5 100.603 1.017 0.338",
        "C C2 PHN 5 99.261 1.077 0.716",
        "C C3 PHN 5 98.403 0.015 0.423",
        "C C4 PHN 5 98.886 -1.108 -0.249",
        "C C5 PHN 5 100.227 -1.170 -0.627",
        "C C6 PHN 5 101.075 -0.108 -0.332",
        "O O1 PHN 5 102.378 -0.209 -0.720",
        "C C1 MAM 6 124.432 0.077 -0.047",
        "N N1 MAM 6 125.856 -0.069 -0.293",
        "C C1 NIT 7 149.514 0.004 -0.008",
        "C C2 NIT 7 150.975 -0.008 0.016",
        "N N1 NIT 7 152.135 -0.018 0.034",
        "C C1 PYR 8 176.143 -0.303 0.058",
        "C C2 PYR 8 175.158 -1.284 -0.005",
        "C C3 PYR 8 173.834 -0.880 -0.072",
        "N N1 PYR 8 173.443 0.413 -0.080",
        "C C4 PYR 8 174.421 1.343 -0.017",
        "C C5 PYR 8 175.771 1.037 0.053",
        "C C1 PRL 9 199.232 -0.905 0.140",
        "C C2 PRL 9 198.899 0.441 -0.144",
        "C C3 PRL 9 200.078 1.148 -0.224",
        "N N1 PRL 9 201.111 0.275 0.002",
        "C C4 PRL 9 200.604 -0.979 0.226",
        "C C1 IMI 10 225.171 0.964 0.224",
        "C C2 IMI 10 226.228 0.123 -0.031",
        "N N1 IMI 10 225.784 -1.142 -0.315",
        "C C3 IMI 10 224.475 -1.077 -0.234",
        "N N2 IMI 10 224.057 0.182 0.091",
        "C C1 MIM 11 248.327 -0.095 -0.015",
        "N N1 MIM 11 249.759 0.049 -0.008",
        "C C2 MIM 11 250.680 -0.954 0.091",
        "C C3 MIM 11 251.899 -0.317 0.054",
        "N N2 MIM 11 251.736 1.040 -0.065",
        "C C4 MIM 11 250.437 1.230 -0.099",
        "C C1 OXZ 12 274.107 -0.467 0.200",
        "C C2 OXZ 12 274.413 0.812 -0.193",
        "O O1 OXZ 12 275.766 0.901 -0.317",
        "C C3 OXZ 12 276.204 -0.346 0.013",
        "N N1 OXZ 12 275.259 -1.202 0.329",
        "C C1 NMA 13 301.863 0.232 0.228",
        "C C2 NMA 13 300.397 0.334 0.550",
        "O O1 NMA 13 299.982 1.000 1.493",
        "N N1 NMA 13 299.582 -0.367 -0.310",
        "C C3 NMA 13 298.148 -0.293 -0.193",
        "O O1 PYO 14 327.610 0.383 -0.060",
        "C C1 PYO 14 326.399 0.195 -0.032",
        "C C2 PYO 14 325.441 1.331 -0.032",
        "C C3 PYO 14 324.123 1.086 -0.001",
        "C C4 PYO 14 323.625 -0.270 0.033",
        "C C5 PYO 14 324.487 -1.295 0.033",
        "N N1 PYO 14 325.838 -1.060 0.001",
        "C C1 MPO 15 348.787 0.239 -0.238",
        "O O1 MPO 15 350.013 -0.167 -0.797",
        "P P1 MPO 15 351.259 -0.262 0.269",
        "O O2 MPO 15 352.430 -0.713 -0.590",
        "O O3 MPO 15 351.375 1.156 0.813",
        "O O4 MPO 15 350.782 -1.289 1.288",
        "C C1 AMD 16 374.285 -0.837 0.028",
        "C C2 AMD 16 375.100 0.355 0.446",
        "N N1 AMD 16 374.559 1.518 0.601",
        "N N2 AMD 16 376.438 0.225 0.633",
        "C C1 IMN 17 398.714 0.347 0.380",
        "C C2 IMN 17 399.346 -0.587 -0.609",
        "N N1 IMN 17 400.602 -0.777 -0.824",
        "C C3 IMN 17 401.569 -0.023 -0.043",
        "C C1 TMA 18 424.701 1.327 0.359",
        "N N1 TMA 18 424.854 -0.114 0.555",
        "C C2 TMA 18 426.254 -0.505 0.398",
        "C C3 TMA 18 424.009 -0.850 -0.384",
        "O O1 CHX 19 451.346 -1.311 -1.000",
        "C C1 CHX 19 451.074 -0.785 0.297",
        "C C2 CHX 19 449.686 -1.250 0.735",
        "C C3 CHX 19 448.574 -0.613 -0.099",
        "C C4 CHX 19 448.681 0.909 -0.110",
        "C C5 CHX 19 450.060 1.374 -0.568",
        "C C6 CHX 19 451.173 0.742 0.267",
        "C C1 SAP 20 476.288 -0.099 0.344",
        "S S1 SAP 20 474.691 -0.510 -0.325",
        "O O1 SAP 20 474.099 -1.569 0.456",
        "O O2 SAP 20 474.769 -0.559 -1.765",
        "N N1 SAP 20 473.806 0.867 0.034",
        "C C1 SAS 21 501.744 -0.144 -0.400",
        "S S1 SAS 21 500.543 0.735 0.578",
        "O O1 SAS 21 500.073 1.883 -0.162",
        "O O2 SAS 21 501.056 0.853 1.924",
        "N N1 SAS 21 499.292 -0.384 0.663",
        "C C2 SAS 21 498.300 -0.380 -0.410",
        "C C1 DMA 22 523.065 0.963 -1.061",
        "N N1 DMA 22 523.687 -0.135 -0.323",
        "C C2 DMA 22 522.826 -0.812 0.644",
        "C C3 DMA 22 525.068 -0.094 -0.088",
        "C C4 DMA 22 525.695 -0.944 0.843",
        "C C5 DMA 22 527.078 -0.918 1.070",
        "C C6 DMA 22 527.886 -0.042 0.359",
        "C C7 DMA 22 527.308 0.800 -0.581",
        "C C8 DMA 22 525.924 0.767 -0.801",
    };
    rigid r = parse_rows(rows);
    struct expectation {
        const char* res;
        int seq;
        const char* name;
        sz ad;
        bool donor;
    };
    const expectation expected[] = {
        {"ETH", 1, "O1", AD_TYPE_OA, true},   // alcohol
        {"ACN", 2, "O1", AD_TYPE_OA, false},  // ketone
        {"ACT", 3, "O1", AD_TYPE_OA, false},  // carboxylate
        {"ACT", 3, "O2", AD_TYPE_OA, false},
        {"ACI", 4, "O1", AD_TYPE_OA, false},  // carboxylic acid: carboxylate at pH 7
        {"ACI", 4, "O2", AD_TYPE_OA, false},
        {"PHN", 5, "O1", AD_TYPE_OA, true},   // phenol
        {"PHN", 5, "C1", AD_TYPE_A, false},   // aromatic ring
        {"MAM", 6, "N1", AD_TYPE_N, true},    // amine
        {"NIT", 7, "N1", AD_TYPE_NA, false},  // nitrile
        {"PYR", 8, "N1", AD_TYPE_NA, false},  // pyridine
        {"PRL", 9, "N1", AD_TYPE_N, true},    // pyrrole
        {"IMI", 10, "N1", AD_TYPE_NA, true},  // imidazole: tautomer unknown
        {"IMI", 10, "N2", AD_TYPE_NA, true},
        {"MIM", 11, "N1", AD_TYPE_N, false},  // N-methyl imidazole N1
        {"MIM", 11, "N2", AD_TYPE_NA, false}, // ... and N3
        {"OXZ", 12, "N1", AD_TYPE_NA, false}, // oxazole
        {"OXZ", 12, "O1", AD_TYPE_OA, false},
        {"NMA", 13, "N1", AD_TYPE_N, true},   // secondary amide
        {"NMA", 13, "O1", AD_TYPE_OA, false},
        {"PYO", 14, "N1", AD_TYPE_N, true},   // 2-pyridone (lactam)
        {"PYO", 14, "O1", AD_TYPE_OA, false},
        {"MPO", 15, "O1", AD_TYPE_OA, false}, // phosphate ester
        {"MPO", 15, "O2", AD_TYPE_OA, false},
        {"MPO", 15, "O3", AD_TYPE_OA, false},
        {"AMD", 16, "N1", AD_TYPE_N, true},   // amidinium
        {"AMD", 16, "N2", AD_TYPE_N, true},
        {"IMN", 17, "N1", AD_TYPE_NA, false}, // imine
        {"TMA", 18, "N1", AD_TYPE_N, true},   // tertiary aliphatic amine: protonated at pH 7
        {"CHX", 19, "O1", AD_TYPE_OA, true},  // cyclohexanol
        {"CHX", 19, "C2", AD_TYPE_C, false},  // non-planar ring: not aromatic
        {"SAP", 20, "N1", AD_TYPE_NA, true},  // primary sulfonamide
        {"SAS", 21, "N1", AD_TYPE_NA, true},  // secondary sulfonamide
        {"DMA", 22, "N1", AD_TYPE_N, false},  // N,N-dimethylaniline: not a basic amine
    };
    for (const expectation& e : expected) {
        const sz i = row_index(rows, e.res, e.seq, e.name);
        INFO(e.res << " " << e.name);
        CHECK(r.atoms[i].ad == e.ad);
        CHECK(r.implicit_donors[i] == e.donor);
    }
}

TEST_CASE("mmcif histidine tautomers from hydrogen-bond partners", "[parse_mmcif_templates]") {
    const std::vector<std::string> rows = {
        "C CB HIS 101 1.893 40.035 0.227",
        "C CG HIS 101 0.422 39.954 0.102",
        "C CD2 HIS 101 -0.426 40.679 -0.707",
        "N NE2 HIS 101 -1.679 40.208 -0.428",
        "C CE1 HIS 101 -1.551 39.234 0.519",
        "N ND1 HIS 101 -0.295 39.057 0.859",
        "O OD1 ASP 102 0.647 37.177 2.782",
        "C CB HIS 103 1.893 80.035 0.227",
        "C CG HIS 103 0.422 79.954 0.102",
        "C CD2 HIS 103 -0.426 80.679 -0.707",
        "N NE2 HIS 103 -1.679 80.208 -0.428",
        "C CE1 HIS 103 -1.551 79.234 0.519",
        "N ND1 HIS 103 -0.295 79.057 0.859",
        "N NZ LYS 104 0.647 77.177 2.782",
        "C CB HIS 105 1.893 120.035 0.227",
        "C CG HIS 105 0.422 119.954 0.102",
        "C CD2 HIS 105 -0.426 120.679 -0.707",
        "N NE2 HIS 105 -1.679 120.208 -0.428",
        "C CE1 HIS 105 -1.551 119.234 0.519",
        "N ND1 HIS 105 -0.295 119.057 0.859",
        "O OD1 ASP 106 0.647 117.177 2.782",
        "O OE1 GLU 107 -4.116 121.095 -1.609",
        "C CB HIS 108 1.893 160.035 0.227",
        "C CG HIS 108 0.422 159.954 0.102",
        "C CD2 HIS 108 -0.426 160.679 -0.707",
        "N NE2 HIS 108 -1.679 160.208 -0.428",
        "C CE1 HIS 108 -1.551 159.234 0.519",
        "N ND1 HIS 108 -0.295 159.057 0.859",
        "O OD1 ASP 109 -0.577 161.015 2.911",
    };
    rigid r = parse_rows(rows);
    const auto check = [&](int seq, bool nd1_h, bool ne2_h) {
        const sz nd1 = row_index(rows, "HIS", seq, "ND1"), ne2 = row_index(rows, "HIS", seq, "NE2");
        INFO("HIS " << seq);
        CHECK(r.implicit_donors[nd1] == nd1_h);
        CHECK(r.implicit_donors[ne2] == ne2_h);
        CHECK(r.atoms[nd1].ad == (nd1_h ? AD_TYPE_N : AD_TYPE_NA));
        CHECK(r.atoms[ne2].ad == (ne2_h ? AD_TYPE_N : AD_TYPE_NA));
    };
    check(101, true, false);  // carboxylate O in front of ND1: HID
    check(103, false, true);  // lysine NZ in front of ND1: HIE
    check(105, true, true);   // carboxylates in front of both: HIP
    check(108, true, false);  // partner not along the N-H direction: default HID
}
