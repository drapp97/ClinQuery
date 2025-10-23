import streamlit as st
from Bio import Entrez
import pandas as pd
import xml.etree.ElementTree as ET

st.title("🧬 ClinQuery")
st.write("Search ClinVar for variant information by gene and HGVS cDNA notation.")

# Set your email for NCBI Entrez
Entrez.email = "daniel.rappaport@sickkids.ca"  # <-- replace with your email

st.sidebar.header("Variant Lookup")

gene = st.sidebar.text_input("Gene symbol (e.g., NPC1):")
variant = st.sidebar.text_input("HGVS cDNA (e.g., c.3044G>T):")

if st.sidebar.button("Search"):
    if not gene or not variant:
        st.warning("Please enter both a gene and a variant.")
    else:
        query = f'"{gene}"[GENE] AND "{variant}"[VARNAME]'
        st.info(f"Searching ClinVar for: **{gene} {variant}** ...")
        try:
            # Step 1: search ClinVar
            handle = Entrez.esearch(db="clinvar", term=query, retmax=5)
            record = Entrez.read(handle)
            handle.close()
            id_list = record["IdList"]

            if not id_list:
                st.warning("No results found.")
            else:
                # Step 2: fetch detailed records
                handle = Entrez.efetch(db="clinvar", id=",".join(id_list), retmode="xml")
                data = handle.read()
                handle.close()

                root = ET.fromstring(data)
                results = []

                for clinvar_record in root.findall(".//ClinVarSet"):
                    significance = clinvar_record.findtext("ReferenceClinVarAssertion/ClinicalSignificance/Description")
                    condition = clinvar_record.findtext("ReferenceClinVarAssertion/TraitSet/Trait/Name/ElementValue")
                    review_status = clinvar_record.findtext("ReferenceClinVarAssertion/ClinicalSignificance/ReviewStatus")
                    rcv_id = clinvar_record.findtext("ClinVarAccession/Acc")
                    results.append({
                        "RCV ID": rcv_id or "N/A",
                        "Clinical Significance": significance or "N/A",
                        "Condition": condition or "N/A",
                        "Review Status": review_status or "N/A"
                    })

                df = pd.DataFrame(results)
                st.success(f"Found {len(df)} matching record(s).")
                st.dataframe(df)
        except Exception as e:
            st.error(f"An error occurred: {e}")
