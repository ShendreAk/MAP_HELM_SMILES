from rdkit import Chem
import re
import copy
import pandas as pd
import warnings

warnings.filterwarnings('ignore')

df_monomers = pd.read_csv('data/MAP_momomers_library_new.csv')

# Copyright (c) 2021-2024 Charles Xu and others
#
# Permission is hereby granted, free of charge, to any person obtaining
# a copy of this software and associated documentation files (the
# "Software"), to deal in the Software without restriction, including
# without limitation the rights to use, copy, modify, merge, publish,
# distribute, sublicense, and/or sell copies of the Software, and to
# permit persons to whom the Software is furnished to do so, subject to
# the following conditions:
#
# The above copyright notice and this permission notice shall be
# included in all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
# EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
# MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
# NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
# LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
# OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
# WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

# code under MIT licence Copyright (c) 2021-2024 Charles Xu and others, starts here with required modifications



def relabel_rgroup2index(smi):
    """
        Input: SMILES of a molecule with R groups, e.g. 'CCCC[C@H](N(C)[*:_R1])C([*:_R2])=O'
        Output: SMILES of a molecule with R groups, e.g. 'CCCC[C@H](N(C)[*:1])C([*:2])=O'
    """
    r_group_name = re.findall(r'\[\*\:(_R\d)\]', smi)
    r_group_name = list(r_group_name)
    # logger.debug(f"r_group_name: {r_group_name}")
    
    checked_rgroup = []
    for name in r_group_name:
        # logger.debug(f"name: {name}, name[2:]: {name[2:]}, checked_rgroup: {checked_rgroup}")
        if name in checked_rgroup:  # For the case of _R3 and _R3, to *:3 and *:4
            # logger.debug(f"name in checked_rgroup: {name}")
            smi = smi.replace(name, f'{int(name[2:])+1}', 1)  # Replace the first occurence
            name = f'_R{int(name[2:])+1}'
        else:  # For the case of _R1 and _R2, to *:1 and *:2
            # logger.debug(f"name not in checked_rgroup: {name}")
            smi = smi.replace(name, f'{name[2:]}', 1)  # Replace the first occurence
            # logger.debug(f"smi: {smi}")
        checked_rgroup.append(name)
        # logger.debug(f"smi: {smi}")
    return smi

def relabel_rgroup2label(smi):
    """
        Input: A SMILES with R groups, e.g. 'CCCC[C@H](N(C)[*:1])C([*:2])=O'
        Output: A SMILES with R groups, e.g. 'CCCC[C@H](N(C)[*:_R1])C([*:_R2])=O'
    """
    r_group_name = re.findall(r'\[(\*\:\d)\]', smi)
    r_group_name = list(r_group_name)
    for name in r_group_name:
        smi = smi.replace(name, f'*:_R{name[2:]}')
    return smi

def get_smi_from_cxsmiles(cxsmiles):
    """
    Get SMILES from CXSMILES
    Input: 'CCCC[C@H](N(C)[*])C([*])=O |$;;;;;;;_R1;;_R2;$|'
    Output: 'CCCC[C@H](N(C)[*:_R1])C([*:_R2])=O'
    """
    smi_list = cxsmiles.split('|')
    smi, pos = smi_list[0], smi_list[1]
    labels = pos.split('$')[1].split(';')

    # Replace * with *:label according to the occurence of * in the SMILES
    for label in labels:
        if len(label) == 0:
            continue
        index = smi.index('[*]')
        smi = smi[:index] + f'[*:{label}]' + smi[index+3:]

    return smi.strip()

def get_cxsmiles_from_smi(smi):
    """
    Get CXSMILES from SMILES
    Input: 'CCCC[C@H](N(C)[*:_R1])C([*:_R2])=O'
    Output: 'CCCC[C@H](N(C)[*])C([*])=O|$;;;;;;;_R1;;_R2;$|'
    """
    cxsmiles = smi
    # Get all labels with pattern [*:label]
    labels = re.findall(r'\[\*\:(.*?)\]', smi)
    r_groups = []
    for label in labels:
        # Replace [*:label] with [*]
        cxsmiles = cxsmiles.replace(f'[*:{label}]', '[*]')
        r_groups.append(f'{label}')     # Record label to r_groups list

    pos = list()
    r_group_idx = 0
    # Iterate through the SMILES and add ';'s in pos
    for i in range(len(cxsmiles)):
        if cxsmiles[i] in ('H', '@', '[', ']', '(', ')', '=', '-', '#', ':', '+',
                           '1', '2', '3', '4', '5', '6', '7', '8', '9', '0', '/', '\\',
                           'l', 'r'):  # 'l' for Cl, 'r' for 'Br'
            # print(cxsmiles[i])
            continue
        elif cxsmiles[i] == '*':
            if r_group_idx >= len(r_groups):
                raise Exception(
                    f'Erorr in converting SMILES to CXSMILES, {smi}')
            pos.append(f"{r_groups[r_group_idx]}")
            r_group_idx += 1
        else:
            pos.append('')
    pos = '|$' + ';'.join(pos) + '$|'
    return f'{cxsmiles} {pos}'

def combine_monomer_unused_rgroup(smi, r_group):
    """ 
    Combine two molecules with unused R groups
    smi: SMILES of the monomer, e.g. 'CCCC[C@H](N(C)[*:1])C([*:2])=O'
    r_group: SMILES of the second molecule, e.g. '[*:1][H]', '[*:2][OH]'

    Note: the R group should be in the form of [*:1], [*:2], etc. Not [*:_R1], [*:_R2], etc.
    """
    # [*]
    m1 = Chem.MolFromSmiles(smi)
    m2 = Chem.MolFromSmiles(r_group)

    mol = Chem.molzip(m1, m2)
    return Chem.MolToSmiles(mol)


def replace_unused_r_groups(mol_smi, r_groups: dict, used_r_groups: list):
    # logger.debug(f"before replace unused r group mol_smi: {mol_smi}, r_groups: {r_groups}, used_r_groups: {used_r_groups}")
    mol_smi = relabel_rgroup2index(mol_smi)
    # logger.debug(f"mol_smi: {mol_smi}")
    for r_group in r_groups.keys():
        if r_group not in used_r_groups:
            # e.g. '[H][*:1]'
            r_group_smi = f"[{r_groups[r_group]}][*:{r_group[1:]}]"
            # logger.debug(f"mol_smi: {mol_smi}, r_group_smi: {r_group_smi}")
            mol_smi = combine_monomer_unused_rgroup(mol_smi, r_group_smi)
    # Remove [H] in smiles, as it is useless
    mol_smi = mol_smi.replace('[H]', '')
    # logger.debug(f"after replace unused r group mol_smi: {mol_smi}")
    return relabel_rgroup2label(mol_smi)

def clean_dummy_labels_in_cxsmiles(smi):
    """
    Clean dummy labels in cxsmiles
    Input: '*C(=O)[C@@H]1CCCN1C(C)=O |$_R2;;;;;;;;;;$,atomProp:0.dummyLabel.*|'
    Output: '[*]C(=O)[C@@H]1CCCN1C(C)=O |$_R2;;;;;;;;;;$|'
    """
    smi_parts = smi.split('|')
    # logger.info(smi_parts)
    smi = smi_parts[0].replace('*', '[*]') + '|$' + \
        smi_parts[1].split('$')[1] + '$|'
    return smi
def combine_fragments(smi1, smi2):
    # logger.info(f"Combine {smi1} and {smi2}...")

    m1 = Chem.MolFromSmiles(get_cxsmiles_from_smi(smi1))
    m2 = Chem.MolFromSmiles(get_cxsmiles_from_smi(smi2))

    for atm in m1.GetAtoms():
        if atm.HasProp("atomLabel") and atm.GetProp("atomLabel") == "_R2":
            atm.SetAtomMapNum(10)  # Use 10 as the atom map number
    for atm in m2.GetAtoms():
        if atm.HasProp("atomLabel") and atm.GetProp("atomLabel") == "_R1":
            atm.SetAtomMapNum(10)  # Use 10 as the atom map number
    mol = Chem.molzip(m1, m2)

    smi = Chem.MolToCXSmiles(mol)
    
    if '|' in smi:
        # logger.info(smi)
        smi = get_smi_from_cxsmiles(clean_dummy_labels_in_cxsmiles(smi))
    # logger.info(f"Combined smiles: {smi}")
    return smi

def get_linear_peptide(monomer_smis):
    """Get linear peptide from monomers"""

    for idx, monomer in enumerate(monomer_smis):
        if idx == 0:
            smi = monomer
        else:
            smi = combine_fragments(smi, monomer)
    return smi

def connect_mapped_atoms(smi, end1, end2):
    ''' 
        smi: smiles contain the two merge ends, e.g. '[*:3]O[C@H](C)[C@H](NC(=O)[C@@H](CC(C)C)N(C)C(=O)[C@@H]1CCCN1C(C)=O)C(=O)N(C)[C@@H](CC(C)C)C(=O)N[C@@H](CC(C)C)C(=O)N[C@@H](CC(C)C)C(=O)N(C)[C@@H](C)C(=O)N[C@H](CC(C)C)C(=O)N1CCC[C@@H]1C([*:2])=O'
        end1: atom map number of the first end, e.g. 3
        end2: atom map number of the second end, e.g. 2

        Note: this will not get atomic stereochemistry right
    '''
    res = Chem.RWMol(Chem.MolFromSmiles(smi))

    dummy1 = None
    dummy2 = None

    for atom in res.GetAtoms():
        if atom.GetAtomMapNum() == end1:
            dummy1 = atom
        elif atom.GetAtomMapNum() == end2:
            dummy2 = atom

    assert dummy1 is not None and dummy2 is not None
    assert dummy1.GetDegree() == 1
    assert dummy2.GetDegree() == 1

    nbr1 = dummy1.GetNeighbors()[0]
    nbr2 = dummy2.GetNeighbors()[0]

    res.BeginBatchEdit()
    res.RemoveAtom(dummy1.GetIdx())
    res.RemoveAtom(dummy2.GetIdx())

    res.AddBond(nbr1.GetIdx(), nbr2.GetIdx(), Chem.BondType.SINGLE)
    res.CommitBatchEdit()
    return Chem.MolToSmiles(res)

def cyclize_linpep_from_smi(smi, link):
    """
    Cyclize a peptide through a link, through definiting a user-defined reaction
    smi: SMILES of a linear peptide, e.g. '[*:_R3]O[C@H](C)[C@H](NC(=O)[C@@H](CC(C)C)N(C)C(=O)[C@@H]1CCCN1C(C)=O)C(=O)N(C)[C@@H](CC(C)C)C(=O)N[C@@H](CC(C)C)C(=O)N[C@@H](CC(C)C)C(=O)N(C)[C@@H](C)C(=O)N[C@H](CC(C)C)C(=O)N1CCC[C@@H]1C([*:_R2])=O'
    link: the link between the two ends, e.g. '4:R3-10:R2', what if '4:R3-14:R3'
    """
    smi = relabel_rgroup2index(smi) # *:_R3 -> *:3
    link = link.split('-')
    start, end = (int(end.split(':')[1][1:]) for end in link) # 3, 3
    if end == start: end = end + 1 # 3, 3 -> 3, 4
    smi = connect_mapped_atoms(smi, start, end)
    return smi

def get_links_between_monomers(monomer_smis, cyclic_link=None):
    monomer_links = {}

    def add_link(source_idx, source_r_group, target_idx, target_r_group):
        # Add link to source monomer
        if source_idx not in monomer_links:
            monomer_links[source_idx] = {source_r_group: (target_idx, target_r_group)}
        else:
            monomer_links[source_idx][source_r_group] = (target_idx, target_r_group)

        # Add link to target monomer
        if target_idx not in monomer_links:
            monomer_links[target_idx] = {target_r_group: None}
        else:
            monomer_links[target_idx][target_r_group] = None

    # Add linear links to monomer_links: R2-R1
    for idx in range(len(monomer_smis)-1):
        add_link(idx, '_R2', idx+1, '_R1')

    # Add cyclic link to monomer_links: '4:R3-10:R2'
    def get_idx_rgroup(node):
        idx, r_group = node.split(':')
        return int(idx) - 1, f'_{r_group}'
    
    if cyclic_link:
        source, target = cyclic_link.split('-')
        source_idx, source_r_group = get_idx_rgroup(source)
        target_idx, target_r_group = get_idx_rgroup(target)
        add_link(source_idx, source_r_group, target_idx, target_r_group)

    return monomer_links

def restore_unused_rgroup(monomer_smis, monomer_r_groups, monomer_links):
    # Replace unused R groups with default _R1, _R2, _R3
    for idx, mol_cxsmi in enumerate(monomer_smis):
        # logger.debug(f"mol_cxsmi: {mol_cxsmi}")
        r_groups = monomer_r_groups[idx]
        used_r_groups = [r_group[1:] for r_group in list(monomer_links[idx].keys())]
        # logger.debug(f"used_r_groups: {used_r_groups}, r_groups: {r_groups}")
        mol_smi = replace_unused_r_groups(mol_cxsmi, r_groups, used_r_groups)
        # logger.debug(f"after replace rgroup mol_smi: {mol_smi}")
        monomer_smis[idx] = mol_smi
    return monomer_smis


monomers2smi_dict = {}
monomers2r_groups_dict = {}
for index, row in df_monomers.iterrows():
    smi = get_smi_from_cxsmiles(row['CXSMILES'])
    monomers2smi_dict[row['Symbol']] = smi
    monomers2r_groups_dict[row['Symbol']] = {}
    for r_group in ['R1', 'R2', 'R3']:
        if row[r_group] != '-':
            monomers2r_groups_dict[row['Symbol']][r_group] = row[r_group]

def cyclize_linpep_from_map(monomer_list, cyclic_link):
    # default monomers with R group denoted by R1, R2, R3
    monomer_smis = [copy.deepcopy(monomers2smi_dict[monomer]) for monomer in monomer_list]
    # Default R groups if no linker
    monomer_r_groups = [copy.deepcopy(monomers2r_groups_dict[monomer]) for monomer in monomer_list]
    
    monomer_links = get_links_between_monomers(monomer_smis, cyclic_link)
    
    monomer_smis = restore_unused_rgroup(monomer_smis, monomer_r_groups, monomer_links)

    # Merge the monomers one by one into a linear peptide
    pep_smi = get_linear_peptide(monomer_smis)
    
    if cyclic_link: # Cyclize the linear peptide
        pep_smi = cyclize_linpep_from_smi(pep_smi, cyclic_link)

    return pep_smi

def linpep_from_map(monomer_list):
    # default monomers with R group denoted by R1, R2, R3
    monomer_smis = [copy.deepcopy(monomers2smi_dict[monomer]) for monomer in monomer_list]
    # Default R groups if no linker
    monomer_r_groups = [copy.deepcopy(monomers2r_groups_dict[monomer]) for monomer in monomer_list]
    
    monomer_links = get_links_between_monomers(monomer_smis, cyclic_link=None)
    
    monomer_smis = restore_unused_rgroup(monomer_smis, monomer_r_groups, monomer_links)

    # Merge the monomers one by one into a linear peptide
    pep_smi = get_linear_peptide(monomer_smis)

    return pep_smi


def extract_data(input_string):
    # Define the regex pattern to match the {cyc: x-y} format with positive integers
    cyc_pattern = r'\{cyc:\s*([N]|\d+)-([C]|\d+)\}'
    
    # Check if the pattern exists in the input string
    match = re.search(cyc_pattern, input_string)
    if match:
        # Store the matched pattern in linker_info
        linker_info = match.group(0)  # This captures the entire matched string
        
        # If found, extract everything after the {cyc: x-y} pattern
        result = re.sub(cyc_pattern, '', input_string).strip()

        nterm_pattern = r'\{nt:[^}]+\}'
        string = ''
        # Find all nterm modifications
        nterm_modifications = re.findall(nterm_pattern, result)
        if nterm_modifications:
            # Remove nterm modifications from final_output
            result = re.sub(nterm_pattern, '', result)
            
            # create string
            string += ''.join(nterm_modifications) + result
            result = f'{string}'
            
            return result, linker_info
        else:
            return result, linker_info
    else:
        nterm_pattern = r'\{nt:[^}]+\}'
        string = ''
        # Find all nterm modifications
        nterm_modifications = re.findall(nterm_pattern, input_string)
        if nterm_modifications:
            # Remove nterm modifications from final_output
            result = re.sub(nterm_pattern, '', input_string)
            
            # create string
            string += ''.join(nterm_modifications) + result
            result = f'{string}'
            
            return result, None
        else:
            return input_string, None  # Return None if the pattern is not found

df2 = pd.read_csv('data/MAP_momomers_library_new.csv')
map_to_helm_dict = df2.set_index('MAP_denotion')['Symbol'].sort_index(ascending=False).to_dict()

def monomer_list_from_linear_seq(linear_seq):
    tokens = []
    i = 0
    while i < len(linear_seq):
        matched = False
        for key in map_to_helm_dict.keys():
            if linear_seq[i:].startswith(key):
                token = map_to_helm_dict[key]
                tokens.append(token)
                i += len(key)
                matched = True
                break
        if not matched:
            # Add unmatched character as-is (can help preserve syntax like {, }, :, etc.)
            tokens.append(linear_seq[i])
            i += 1
    return tokens

def get_smi_from_map(map):
    linear_seq, linker = extract_data(map)
    print(linear_seq)
    # print(linker)
    monomer_list = monomer_list_from_linear_seq(linear_seq)
    if linker:
        print(linker)
        linker_list = linker.split('-')
        cyclic_linker = ''
        if linker_list[0][-1] == '1' :
            cyclic_linker += f'{1}:R1-'
        elif linker_list[0][-1] == 'N' :
            cyclic_linker += f'{1}:R1-'
        else:
            cyclic_linker += f'{int(linker_list[0][-1])}:R3-'
        end_conn = linker_list[1].replace('}','')
        if end_conn == 'C':
            cyclic_linker += f'{len(monomer_list)}:R2'
        elif  end_conn == str(len(monomer_list)):
            cyclic_linker += f'{len(monomer_list)}:R2'
        else:
            cyclic_linker += f'{int(end_conn)}:R3'
        try:
            print(cyclic_linker)
            smi = cyclize_linpep_from_map(monomer_list, cyclic_linker)
            return smi
        except Exception as e:
            return None
    else:
        try:
            smi = linpep_from_map(monomer_list)
            return smi
        except Exception as e:
            return None
          
# code under MIT licence Copyright (c) 2021-2024 Charles Xu and others, ends here 

##HELM to MAP
# Load the MAP monomer library
df2 = pd.read_csv('data/MAP_momomers_library_new.csv')
map_to_helm_dict = df2.set_index('MAP_denotion')['Symbol'].sort_index(ascending=False).to_dict()

def helm_to_map(helm):
    try:
        start = helm.index('{') + 1
        end = helm.index('}')
        helm_sequence = helm[start:end]
        elements = [elem.strip('[]') for elem in helm_sequence.split('.')]
        num_elements = len(elements)
        map_format = ''
        for element in elements:
            map_notation = df2.loc[df2['Symbol'] == element, 'MAP_denotion']
            if not map_notation.empty:
                map_format += map_notation.values[0]
        dollar_split = helm.split('$')
        if len(dollar_split) > 2:
            last_part = dollar_split[1]
            if len(last_part) > 0:
                last_element = last_part.split(',')[-1]
                first_part = last_element.split(':')[0]
                second_part = last_element.split(':')[1].split('-')[1]

                if int(first_part) == 1 and int(second_part)==num_elements:
                    cyc_string = f'N-C'
                else:
                    cyc_string = f'{first_part}-{second_part}'
                final_output = f'{map_format}'
                nterm_pattern = r'\{nt:[^}]+\}'
                cterm_pattern = r'\{ct:[^}]+\}'
                nterm_modifications = re.findall(nterm_pattern, final_output)
                final_output = re.sub(nterm_pattern, '', final_output)
                cterm_modifications = re.findall(cterm_pattern, final_output)
                final_output = re.sub(cterm_pattern, '', final_output)
                final_output += ''.join(nterm_modifications) + ''.join(cterm_modifications)
                return f'{final_output}{{cyc:{cyc_string}}}'
            else:
                final_output = f'{map_format}'
                nterm_pattern = r'\{nt:[^}]+\}'
                cterm_pattern = r'\{ct:[^}]+\}'
                nterm_modifications = re.findall(nterm_pattern, final_output)
                final_output = re.sub(nterm_pattern, '', final_output)
                cterm_modifications = re.findall(cterm_pattern, final_output)
                final_output = re.sub(cterm_pattern, '', final_output)
                final_output += ''.join(nterm_modifications) + ''.join(cterm_modifications)
                return final_output
    except Exception as e:
        return f"ERROR: {e}"
    return ''


##MAP to HELM sequence
def process_HELM_seq(helm_seq, ID):
    if '{' in helm_seq:
        start = helm_seq.index('{')
        end = helm_seq.index('}')
        cyc_seq = helm_seq[start:end]
        seq_len = len(helm_seq[end+1:].split('.'))
        cyc_list = cyc_seq.split('-')
        start_pos = cyc_list[0][-1]
        end_pos = cyc_list[1]
       
        if start_pos == 'N' and end_pos == 'C':
            return f'PEPTIDE{ID}{{{helm_seq[end+1:]}}}$PEPTIDE{ID},PEPTIDE{ID},1:R1-{seq_len}:R2$$$'
        elif start_pos != '1' and end_pos == str(seq_len):
            return f'PEPTIDE{ID}{{{helm_seq[end+1:]}}}$PEPTIDE{ID},PEPTIDE{ID},{start_pos}:R3-{seq_len}:R2$$$'
        elif start_pos == '1' and end_pos != str(seq_len):
            return f'PEPTIDE{ID}{{{helm_seq[end+1:]}}}$PEPTIDE{ID},PEPTIDE{ID},1:R1-{end_pos}:R3$$$'
        else:
            return f'PEPTIDE{ID}{{{helm_seq[end+1:]}}}$PEPTIDE{ID},PEPTIDE{ID},{start_pos}:R3-{end_pos}:R3$$$'
    else:
        return f'PEPTIDE{ID}{{{helm_seq}}}$$$$'

def convert_map_to_helm_sequence(map_str, ID):
    nterm_pattern = r'\{nt:[^}]+\}'
    cyc_pattern = r'\{cyc:\s*([N]|\d+)-([C]|\d+)\}'
    string = ''
    nterm_modifications = re.findall(nterm_pattern, map_str)
    map_str = re.sub(nterm_pattern, '', map_str)
    cyc_string = re.search(cyc_pattern, map_str)
    # print("cyc_string",cyc_string[0])
    map_str = re.sub(cyc_pattern, '', map_str)
    if cyc_string:
        string += ''.join(cyc_string[0]) + ''.join(nterm_modifications) + map_str
    else:
        string += ''.join(nterm_modifications) + map_str

    tokens = []
    i = 0
    while i < len(string):
        matched = False
        for key in map_to_helm_dict.keys():
            if string[i:].startswith(key):
                if string[i-4:i] == 'cyc:':
                    val = map_to_helm_dict[key]
                    token = f'[{val}]' if len(val) > 1 else f'{val}'
                    tokens.append(token)
                    i += len(key)
                    matched = True
                    break
                else:
                    val = map_to_helm_dict[key]
                    token = f'[{val}].' if len(val) > 1 else f'{val}.'
                    tokens.append(token)
                    i += len(key)
                    matched = True
                    break

        if not matched:
            tokens.append(string[i])
            i += 1
    helm_seq = ''.join(tokens).rstrip('.')
    # print(helm_seq)
    return helm_seq
