import argparse
import os
import pandas as pd
import re
from utils import get_smi_from_map, helm_to_map, process_HELM_seq, convert_map_to_helm_sequence 





def main():
    parser = argparse.ArgumentParser(description="HELM-MAP-SMILES Format Converter")
    parser.add_argument("mode", choices=["helm_to_map", "map_to_helm", "map_to_smiles"], help="Conversion mode")
    parser.add_argument("--input", required=True, help="Input string or input file path")
    parser.add_argument("--output", help="Output file path (required if input is a file)")
    parser.add_argument("--id", help="Peptide ID(s) for MAP to HELM (comma-separated for batch)")

    args = parser.parse_args()
    is_file = os.path.isfile(args.input)

    # Handle HELM to MAP
    if args.mode == "helm_to_map":
        if is_file:
            with open(args.input, "r") as f:
                lines = [line.strip() for line in f if line.strip()]
            results = [helm_to_map(line) for line in lines]
            if not args.output:
                raise ValueError("Output file path required for file input.")
            with open(args.output, "w") as f:
                f.write("\n".join(results))
            print(f"Conversion complete. Output saved to {args.output}")
        else:
            print(helm_to_map(args.input))

    # Handle MAP to HELM
    elif args.mode == 'map_to_helm':
      if is_file:
          if not args.output:
              print("Error: Please specify --output for saving the converted HELM format.")
          else:
              with open(args.input, 'r') as infile, open(args.output, 'w') as outfile:
                  for line in infile:
                      line = line.strip()
                      if not line or ',' not in line:
                          continue
                      try:
                          map_seq, peptide_id = map(str.strip, line.split(',', 1))
                          helm_seq = convert_map_to_helm_sequence(map_seq, peptide_id)
                          result = process_HELM_seq(helm_seq, peptide_id)
                          outfile.write(result + '\n')
                      except Exception as e:
                          print(f"Skipping invalid line: {line} | Error: {e}")
              print(f"Conversion complete. Output saved to {args.output}")
      else:
          if not args.id:
              print("Error: --id is required for single sequence MAP to HELM conversion.")
          else:
              helm_seq = convert_map_to_helm_sequence(args.input, args.id)
              result = process_HELM_seq(helm_seq, args.id)
              print(result)



    # Handle MAP to SMILES
    elif args.mode == "map_to_smiles":
        if is_file:
            with open(args.input, "r") as f:
                lines = [line.strip() for line in f if line.strip()]
                # print("lines",lines)
            results = [get_smi_from_map(line) for line in lines]
            # print(results)
            if not args.output:
                raise ValueError("Output file path required for file input.")
            with open(args.output, "w") as f:
                f.write("\n".join(results))
            print(f"Conversion complete. Output saved to {args.output}")
        else:
            print(get_smi_from_map(args.input))

if __name__ == "__main__":
    main()