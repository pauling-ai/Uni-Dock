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

*/

#ifndef VINA_PARSE_MMCIF_H
#define VINA_PARSE_MMCIF_H

#include <string>
#include "parse_pdbqt.h"

// true for *.cif and *.mmcif (case-insensitive)
bool is_mmcif_file_name(const std::string& name);

// Reads the _atom_site category of the first data block directly into receptor atoms, assigning
// AutoDock types from element, connectivity and (for residues without hydrogens) residue
// templates. Fills r.implicit_donors. Can throw struct_parse_error.
void parse_mmcif_rigid(const path& name, rigid& r);
void parse_mmcif_rigid(std::istream& in, rigid& r);

#endif
