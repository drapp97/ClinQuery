import streamlit as st
from Bio import Entrez
import pandas as pd
import xml.etree.ElementTree as ET

st.title("🧬 ClinQuery")
st.write("Search ClinVar variants by gene. Optional HGVS cDNA filter.")

# Set your email for NCBI Entrez
Entrez.email = "your_email@example.com"  # <-- replace with your real email

st.sidebar.header("Variant Lookup")
gene = st.sidebar.text_input("Gene symbol (e.g., BRCA1):")
variant = st.sidebar.text_input("Optional: HGVS cDNA (e.g., c.68_69delAG)")

if st.sidebar.button("Search"):
    if not gene:
        st.warning("Please enter a gene.")
    else:
        st.info(f"Searching ClinVar for gene: {gene} ...")
        try:
            # Step 1: Search ClinVar by gene to get IDs
            search_handle = Entrez.esearch(db="clinvar", term=f"{gene}[GENE]", retmax=100)
            search_record = Entrez.read(search_handle)
            search_handle.close()
            id_list = search_record["IdList"]

            if not id_list:
                st.warning(f"No variants found for gene {gene}.")
            else:
                # Step 2: Fetch detailed records
                fetch_handle = Entrez.efetch(db="clinvar", id=",".join(id_list), retmode="xml")
                fetch_data = fetch_handle.read()
                fetch_handle.close()

                root = ET.fromstring(fetch_data)
                results = []

                for clinvar_record in root.findall(".//ClinVarSet"):
                    # Safe parsing of core fields
                    rcv_id = clinvar_record.findtext("ClinVarAccession/Acc") or "N/A"
                    clinical_significance = clinvar_record.findtext("ReferenceClinVarAssertion/ClinicalSignificance/Description") or "N/A"
                    review_status = clinvar_record.findtext("ReferenceClinVarAssertion/ClinicalSignificance/ReviewStatus") or "N/A"
                    condition = clinvar_record.findtext("ReferenceClinVarAssertion/TraitSet/Trait/Name/ElementValue") or "N/A"

                    # Collect HGVS strings if present
                    hgvs_c = []
                    hgvs_p = []
                    measures = clinvar_record.findall(".//Measure")
                    for m in measures:
                        for attr in m.findall(".//AttributeSet/Attribute"):
                            t_type = attr.findtext("Type")
                            t_value = attr.findtext("Value")
                            if t_type and t_value:
                                if t_type == "HGVS cDNA":
                                    hgvs_c.append(t_value)
                                elif t_type == "HGVS protein":
                                    hgvs_p.append(t_value)

                    results.append({
                        "RCV ID": rcv_id,
                        "Clinical Significance": clinical_significance,
                        "Review Status": review_status,
                        "Condition": condition,
                        "HGVS cDNA": ", ".join(hgvs_c) if hgvs_c else "N/A",
                        "HGVS Protein": ", ".join(hgvs_p) if hgvs_p else "N/A"
                    })

                df = pd.DataFrame(results)

                # Optional: local HGVS filter
                if variant:
                    df_filtered = df[df["HGVS cDNA"].str.contains(variant, na=False)]
                    if not df_filtered.empty:
                        st.success(f"Found {len(df_filtered)} matching record(s) for HGVS variant {variant}.")
                        st.dataframe(df_filtered)
                    else:
                        st.warning(f"No records exactly match HGVS variant '{variant}'. Showing all variants for gene {gene}:")
                        st.dataframe(df)
                else:
                    st.success(f"Found {len(df)} variants for gene {gene}.")
                    st.dataframe(df)

        except Exception as e:
            st.error(f"An error occurred: {e}")
