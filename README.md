# HELM-MAP-SMILES Format Converter CLI

A command-line tool for converting between:

- **HELM notation** → **MAP format**
- **MAP format** → **HELM notation**
- **MAP format** → **SMILES representation**

---

## 📦 Installation & Requirements

### 1. Create a virtual environment (recommended):

```bash
python -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate
```

### 2. transfer the files and folder to the virtual environment

### 3. Then install dependencies

```bash
pip install -r requirements.txt
```

### 4. Required Files

Ensure the following file is present in the directory:

- `data/MAP_momomers_library_new.csv` – Monomer library mapping MAP to HELM
- `utils.py`

---

## 🚀 Usage

The script is executed via command-line and accepts:

- A single sequence via `--input`
- A file with multiple sequences (`--input`) and optional `--output`

## 🔁 Supported Modes

| Mode            | Description                          | Notes                                       |
| --------------- | ------------------------------------ | ------------------------------------------- |
| `helm_to_map`   | Converts a HELM sequence to MAP      | Input: HELM format                          |
| `map_to_helm`   | Converts MAP format to HELM          | Requires `--id` for peptide IDs             |
| `map_to_smiles` | Converts MAP format to SMILES string | Uses external function `get_smi_from_map()` |

---

```bash
python converter_cli.py <mode> --input <sequence_or_filepath> [--output <output_file>] [--id <peptide_id>]
```

---



## 🧪 Examples

### 🧬 Single Sequence Conversion

```bash
python converter_cli.py helm_to_map --input "PEPTIDE1{[Abu].[Sar].[meL].V.[meL].A.[dA].[meL].[meL].[meV].[Me_Bmt(E)]}$PEPTIDE1,PEPTIDE1,1:R1-11:R2$$$"

python converter_cli.py map_to_helm --input "LL{d}L{d}L{d}PY{cyc:1-6}" --id 1

python converter_cli.py map_to_smiles --input "LL{d}L{d}L{d}PY{cyc:1-6}"
```

---

### 📁 Batch File Conversion

**Input files should have one sequence per line.**

```bash
# HELM to MAP
python converter_cli.py helm_to_map --input input_helm.txt --output output_map.txt

📌 Make sure your input file follows this format:
"PEPTIDE1{[Abu].[Sar].[meL].V.[meL].A.[dA].[meL].[meL].[meV].[Me_Bmt(E)]}$PEPTIDE1,PEPTIDE1,1:R1-11:R2$$$"
"PEPTIDE2{[dL].[dL].L.[dL].P.Y}$PEPTIDE2,PEPTIDE2,1:R1-6:R2$$$"
"PEPTIDE3{[dL].[dL].[dL].[dL].P.Y}$PEPTIDE3,PEPTIDE3,1:R1-6:R2$$$"

# MAP to HELM (IDs required, match the number of lines)
python converter_cli.py --mode map_to_helm --input input_map_ids.txt --output output_helm.txt

📌 Make sure your input file follows this format:
L{d}L{d}LL{d}PY{cyc:1-6},1
LLLL{d}PY{cyc:1-6},2

# MAP to SMILES
python converter_cli.py map_to_smiles --input input_map.txt --output output_smiles.txt

📌 Make sure your input file follows this format:
L{d}L{d}LL{d}PY{cyc:1-6}
L{d}L{d}LL{d}PY{cyc:1-6}
```

---

## 🧱 File Structure

```
.
├── converter_cli.py
├── utils.py
├── requirements.txt
├── README.md
├── app.py
└── data/
    └── MAP_momomers_library_new.csv
```

---

## 🛠 Dependencies

```
rdkit
pandas
streamlit
```

Install them using:

```bash
pip install -r requirements.txt
```

---

# APP
# HELM / MAP / SMILES Format Converter (Streamlit App)

This **Streamlit web app** provides a simple interface to convert peptide sequence formats between:

- **HELM → MAP**
- **MAP → HELM**
- **MAP → SMILES**

---

## 🚀 Features

- Paste or upload sequences in HELM, MAP, or SMILES format
- View conversion results instantly
- Download converted output files
- Includes example sequences for easy testing

---

## 🧠 How to Run the App

Run the Streamlit application using the following command:

```bash
streamlit run app.py
```

This will open the app in your browser.

---

## 📝 Input Format Examples

### 1️⃣ HELM to MAP

**Input:**

```
PEPTIDE1{[Abu].[Sar].[meL].V.[meL].A.[dA].[meL].[meL].[meV].[Me_Bmt(E)]}$PEPTIDE1,PEPTIDE1,1:R1-11:R2$$$
```

---

### 2️⃣ MAP to HELM

**MAP Input:**

```
{nnr:ABU}G{nnm:NMX}L{nnm:NMX}VL{nnm:NMX}AA{d}L{nnm:NMX}L{nnm:NMX}V{nnm:NMX}{nnr:MBM}{cyc:N-C}

```

**Peptide ID:**

```
1
```

---

### 3️⃣ MAP to SMILES

**MAP Input:**

```
L{d}L{d}LL{d}PY{cyc:1-6}
```

---
---

## 🧪 Example Usage

- Navigate to a tab (e.g., "HELM to MAP")
- Paste or load example sequences
- Click "Convert"
- Download the converted output if needed

---


## 📬 Contact

For issues or questions, please open an issue or contact the maintainer.
© 2025 HELM-MAP-SMILES Converter Team
