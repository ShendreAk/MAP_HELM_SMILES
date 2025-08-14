import streamlit as st
import pandas as pd
import re
from utils import get_smi_from_map, helm_to_map, process_HELM_seq, convert_map_to_helm_sequence



st.set_page_config(layout="wide")

def main():
    st.title("MAP/HELM/SMILES Format Converter")

    st.markdown("""
    <style>
    textarea {
        min-height: 300px !important;  /* Taller text areas */
        width: 100% !important;        /* Full width */
    }
    </style>
    """, unsafe_allow_html=True)

    # Tabs for different functionalities
    tab1, tab2, tab3 = st.tabs(["HELM to MAP", "MAP to HELM", "MAP to SMILES"])

    with tab1:
        st.header("Convert HELM to MAP")
        helm_examples = "\n".join([
            "PEPTIDE1{[Abu].[Sar].[meL].V.[meL].A.[dA].[meL].[meL].[meV].[Me_Bmt(E)]}$PEPTIDE1,PEPTIDE1,1:R1-11:R2$$$",
            "PEPTIDE2{[dL].[dL].L.[dL].P.Y}$PEPTIDE2,PEPTIDE2,1:R1-6:R2$$$",
            "PEPTIDE3{[dL].[dL].[dL].[dL].P.Y}$PEPTIDE3,PEPTIDE3,1:R1-6:R2$$$",
            "PEPTIDE4{L.L.L.[dL].P.Y}$PEPTIDE4,PEPTIDE4,1:R1-6:R2$$$",
            "PEPTIDE5{L.[dL].[dL].[dL].P.Y}$PEPTIDE5,PEPTIDE5,1:R1-6:R2$$$",
            "PEPTIDE6{L.L.L.L.[dP].Y}$PEPTIDE6,PEPTIDE6,1:R1-6:R2$$$",
            "PEPTIDE7{[dL].[dL].[dL].[dL].[dP].Y}$PEPTIDE7,PEPTIDE7,1:R1-6:R2$$$",
            "PEPTIDE8{L.L.[dL].[dL].P.Y}$PEPTIDE8,PEPTIDE8,1:R1-6:R2$$$",
            "PEPTIDE9{L.[dL].L.[dL].[dP].Y}$PEPTIDE9,PEPTIDE9,1:R1-6:R2$$$",
            "PEPTIDE10{L.[dL].L.L.[dP].Y}$PEPTIDE10,PEPTIDE10,1:R1-6:R2$$$",
            "PEPTIDE11{[Tyr(Me)].P.[dL].[dL].[dL].L}$PEPTIDE11,PEPTIDE11,1:R1-6:R2$$$",
            "PEPTIDE12{[Tyr(Me)].P.[dP].[dL].L.L}$PEPTIDE12,PEPTIDE12,1:R1-6:R2$$$",
            "PEPTIDE13{[Tyr(Me)].P.L.L.[dP].[dL]}$PEPTIDE13,PEPTIDE13,1:R1-6:R2$$$",
            "PEPTIDE14{[Tyr(Me)].P.L.[dL].[dL].[dL]}$PEPTIDE14,PEPTIDE14,1:R1-6:R2$$$",
            "PEPTIDE15{[Tyr(Me)].P.[dL].[dL].[dL].[dL]}$PEPTIDE15,PEPTIDE15,1:R1-6:R2$$$",
            "PEPTIDE16{[Tyr(Me)].P.[dL].L.[dL].[dL]}$PEPTIDE16,PEPTIDE16,1:R1-6:R2$$$",
            "PEPTIDE17{[Tyr(Me)].[dP].L.[dP].[dL].P.[dL]}$PEPTIDE17,PEPTIDE17,1:R1-7:R2$$$",
            "PEPTIDE18{[Tyr(Me)].[dP].L.L.[dL].L}$PEPTIDE18,PEPTIDE18,1:R1-6:R2$$$",
            "PEPTIDE19{[Tyr(Me)].[dP].L.[dL].P.L.L}$PEPTIDE19,PEPTIDE19,1:R1-7:R2$$$",
            "PEPTIDE20{[Tyr(Me)].[dP].L.[dP].[dL].[dL].L}$PEPTIDE20,PEPTIDE20,1:R1-7:R2$$$",
            "PEPTIDE21{[Tyr(Me)].[dP].L.L.[dL].L.[dL]}$PEPTIDE21,PEPTIDE21,1:R1-7:R2$$$",
            "PEPTIDE48{A.A.L.[meV].L.F.F.P.I.T.G.D.[-pip]}$PEPTIDE48,PEPTIDE48,1:R1-12:R3$$$",
            "PEPTIDE49{A.V.[meA].L.L.I.F.L.P.T.L.D.[-pip]}$PEPTIDE49,PEPTIDE49,1:R1-12:R3$$$",
            "PEPTIDE50{A.I.P.[meL].T.G.F.[meA].V.L.L.D.[-pip]}$PEPTIDE50,PEPTIDE50,1:R1-12:R3$$$",
            "PEPTIDE51{A.L.[Sar].L.[meI].[meA].F.T.F.[meL].P.D.[-pip]}$PEPTIDE51,PEPTIDE51,1:R1-12:R3$$$",
            "PEPTIDE52{A.I.P.[meL].[Sar].T.[meF].[meA].L.[meV].F.D.[-pip]}$PEPTIDE52,PEPTIDE52,1:R1-12:R3$$$",
            "PEPTIDE53{A.G.[meL].[meI].[meA].[meF].T.[meF].[meL].P.D.[-pip]}$PEPTIDE53,PEPTIDE53,1:R1-11:R3$$$",
            "PEPTIDE54{A.L.[meL].F.[meF].I.P.T.L.D.[-pip]}$PEPTIDE54,PEPTIDE54,1:R1-10:R3$$$",
            "PEPTIDE55{A.L.[meI].[meA].F.T.L.[meF].[meV].D.[-pip]}$PEPTIDE55,PEPTIDE55,1:R1-10:R3$$$",
            "PEPTIDE56{A.A.[meL].L.[meI].F.F.P.T.[meL].D.[-pip]}$PEPTIDE56,PEPTIDE56,1:R1-11:R3$$$",
            "PEPTIDE57{A.L.[meI].L.[Sar].T.F.A.[meL].F.V.D.[-pip]}$PEPTIDE57,PEPTIDE57,1:R1-12:R3$$$",
            "PEPTIDE58{A.[meA].[meF].[meL].T.[Sar].[meL].[Ser(tBu)].[meI].D.[-pip]}$PEPTIDE58,PEPTIDE58,1:R1-10:R3$$$",
            "PEPTIDE59{A.[meF].[meL].T.[Sar].[meL].[Ser(tBu)].[meI].D.[-pip]}$PEPTIDE59,PEPTIDE59,1:R1-9:R3$$$",
            "PEPTIDE60{A.[meA].[meF].[meL].T.[Sar].[meL].[Ser(tBu)].D.[-pip]}$PEPTIDE60,PEPTIDE60,1:R1-9:R3$$$"

        ])
        if st.button("Example sequence", key="helm_example_button"):
            st.session_state.helm_input = helm_examples

        helm_input = st.text_area("Enter HELM notations (one per line):", value=st.session_state.get("helm_input", ""), key="helm_input_text")
        if st.button("Convert", key="helm_to_map_convert"):
            if helm_input:
                helm_lines = [line.strip() for line in helm_input.strip().split('\n') if line.strip()]
                map_outputs = [helm_to_map(line) for line in helm_lines]
                map_output = "\n".join(map_outputs)
                st.success("Conversion Successful!")
                st.text_area("MAP Format:", value=map_output,  key="helm_map_output")
                st.download_button(
                    label="Download MAP Format",
                    data=map_output,
                    file_name='output_map.txt',
                    mime='text/plain',
                    key="download_map"
                )
            else:
                st.error("Please enter a valid HELM notation.")

    with tab2:
        st.header("Convert MAP to HELM")
        map_examples = "\n".join([
            "{nnr:ABU}G{nnm:NMX}L{nnm:NMX}VL{nnm:NMX}AA{d}L{nnm:NMX}L{nnm:NMX}V{nnm:NMX}{nnr:MBM}{cyc:N-C}",
            "L{d}L{d}LL{d}PY{cyc:N-C}",
            "L{d}L{d}L{d}L{d}PY{cyc:N-C}",
            "LLLL{d}PY{cyc:N-C}",
            "LL{d}L{d}L{d}PY{cyc:N-C}",
            "LLLLP{d}Y{cyc:N-C}",
            "L{d}L{d}L{d}L{d}P{d}Y{cyc:N-C}",
            "LLL{d}L{d}PY{cyc:N-C}",
            "LL{d}LL{d}P{d}Y{cyc:N-C}",
            "LL{d}LLP{d}Y{cyc:N-C}",
            "Y{nnm:OME}PL{d}L{d}L{d}L{cyc:N-C}",
            "Y{nnm:OME}PP{d}L{d}LL{cyc:N-C}",
            "Y{nnm:OME}PLLP{d}L{d}{cyc:N-C}",
            "Y{nnm:OME}PLL{d}L{d}L{d}{cyc:N-C}",
            "Y{nnm:OME}PL{d}L{d}L{d}L{d}{cyc:N-C}",
            "Y{nnm:OME}PL{d}LL{d}L{d}{cyc:N-C}",
            "Y{nnm:OME}P{d}LP{d}L{d}PL{d}{cyc:N-C}",
            "Y{nnm:OME}P{d}LLL{d}L{cyc:N-C}",
            "Y{nnm:OME}P{d}LL{d}PLL{cyc:N-C}",
            "Y{nnm:OME}P{d}LP{d}L{d}L{d}L{cyc:N-C}",
            "Y{nnm:OME}P{d}LLL{d}LL{d}{cyc:N-C}",
            "AALV{nnm:NMX}LFFPITGD{ct:PPD}{cyc:1-12}",
            "AVA{nnm:NMX}LLIFLPTLD{ct:PPD}{cyc:1-12}",
            "AIPL{nnm:NMX}TGFA{nnm:NMX}VLLD{ct:PPD}{cyc:1-12}",
            "ALG{nnm:NMX}LI{nnm:NMX}A{nnm:NMX}FTFL{nnm:NMX}PD{ct:PPD}{cyc:1-12}",
            "AIPL{nnm:NMX}G{nnm:NMX}TF{nnm:NMX}A{nnm:NMX}LV{nnm:NMX}FD{ct:PPD}{cyc:1-12}"


        ])
        id_examples = "\n".join(["1", "2", "3", "4", "5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26"])
        if st.button("Show Examples", key="map_example_button"):
            st.session_state.map_input = map_examples
            st.session_state.id_input = id_examples

        map_input = st.text_area("Enter MAP sequences (one per line):", value=st.session_state.get("map_input", ""), key="map_input_text")
        id_input = st.text_area("Enter a base ID for the peptide(s):", value=st.session_state.get("id_input", ""), key="map_id_input_text")

        if st.button("Convert", key="map_to_helm_convert"):
            if map_input and id_input:
                map_lines = [line.strip() for line in map_input.strip().split('\n') if line.strip()]
                id_lines = [id.strip() for id in id_input.strip().split('\n') if id.strip()]
                helm_outputs = []
                for line, i in zip(map_lines, id_lines):
                    helm_seq = convert_map_to_helm_sequence(line, i)
                    helm_outputs.append(process_HELM_seq(helm_seq, i))
                helm_output = "\n".join(helm_outputs)
                st.success("Conversion Successful!")
                st.text_area("HELM Format:", value=helm_output,  key="map_helm_output")
                st.download_button(
                    label="Download HELM Format",
                    data=helm_output,
                    file_name='output_helm.txt',
                    mime='text/plain',
                    key="download_helm"
                )
            else:
                st.error("Please enter a valid MAP format and ID.")

    with tab3:
        st.header("Convert MAP to SMILES")
        map_examples = "\n".join([
            "{nnr:ABU}G{nnm:NMX}L{nnm:NMX}VL{nnm:NMX}AA{d}L{nnm:NMX}L{nnm:NMX}V{nnm:NMX}{nnr:MBM}{cyc:N-C}",
            "L{d}L{d}LL{d}PY{cyc:N-C}",
            "L{d}L{d}L{d}L{d}PY{cyc:N-C}",
            "LLLL{d}PY{cyc:N-C}",
            "LL{d}L{d}L{d}PY{cyc:N-C}",
            "LLLLP{d}Y{cyc:N-C}",
            "L{d}L{d}L{d}L{d}P{d}Y{cyc:N-C}",
            "LLL{d}L{d}PY{cyc:N-C}",
            "LL{d}LL{d}P{d}Y{cyc:N-C}",
            "LL{d}LLP{d}Y{cyc:N-C}",
            "Y{nnm:OME}PL{d}L{d}L{d}L{cyc:N-C}",
            "Y{nnm:OME}PP{d}L{d}LL{cyc:N-C}",
            "Y{nnm:OME}PLLP{d}L{d}{cyc:N-C}",
            "Y{nnm:OME}PLL{d}L{d}L{d}{cyc:N-C}",
            "Y{nnm:OME}PL{d}L{d}L{d}L{d}{cyc:N-C}",
            "Y{nnm:OME}PL{d}LL{d}L{d}{cyc:N-C}",
            "Y{nnm:OME}P{d}LP{d}L{d}PL{d}{cyc:N-C}",
            "Y{nnm:OME}P{d}LLL{d}L{cyc:N-C}",
            "Y{nnm:OME}P{d}LL{d}PLL{cyc:N-C}",
            "Y{nnm:OME}P{d}LP{d}L{d}L{d}L{cyc:N-C}",
            "Y{nnm:OME}P{d}LLL{d}LL{d}{cyc:N-C}",
            "AALV{nnm:NMX}LFFPITGD{ct:PPD}{cyc:1-12}",
            "AVA{nnm:NMX}LLIFLPTLD{ct:PPD}{cyc:1-12}",
            "AIPL{nnm:NMX}TGFA{nnm:NMX}VLLD{ct:PPD}{cyc:1-12}",
            "ALG{nnm:NMX}LI{nnm:NMX}A{nnm:NMX}FTFL{nnm:NMX}PD{ct:PPD}{cyc:1-12}",
            "AIPL{nnm:NMX}G{nnm:NMX}TF{nnm:NMX}A{nnm:NMX}LV{nnm:NMX}FD{ct:PPD}{cyc:1-12}"
        ])

        if st.button("Example sequence", key="smiles_example_button"):
            st.session_state.map_input_smiles = map_examples

        map_input_smiles = st.text_area("Enter MAP sequences (one per line):", value=st.session_state.get("map_input_smiles", ""), key="map_smiles_input_text")

        if st.button("Convert", key="map_to_smiles_convert"):
            if map_input_smiles:
                map_lines = [line.strip() for line in map_input_smiles.strip().split('\n') if line.strip()]
                smiles_outputs = [get_smi_from_map(line) for line in map_lines]
                print(smiles_outputs)
                smiles_output = "\n".join(smiles_outputs)
                st.success("Conversion Successful!")
                st.text_area("SMILES Format:", value=smiles_output,  key="map_smiles_output")
                st.download_button(
                    label="Download SMILES Format",
                    data=smiles_output,
                    file_name='output_smiles.txt',
                    mime='text/plain',
                    key="download_smiles"
                )
            else:
                st.error("Please enter a valid MAP format.")


        

if __name__ == "__main__":
    main()