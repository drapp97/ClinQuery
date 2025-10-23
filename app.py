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
    if not gene:
        st.warning("Please enter a gene.")
    else:
        st.info(f"Searching ClinVar for gene: **{gene}** ...")
        try:
            # Step 1: search ClinVar by gene only
            query = f'"{gene}"[GENE]'
            handle = Entrez.esearch(db="clinvar", term=query, retmax=50)  # fetch up to 50 variants
            record = Entrez.read(handle)
            handle.close()
            id_list = record["IdList"]

            if not id_list:
                st.warning("No results found for that gene.")
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

                    # Get HGVS strings (cDNA / protein)
                    hgvs_c = []
                    hgvs_p = []
                    measures = clinvar_record.findall(".//Measure")
                    for m in measures:
                        for t in m.findall(".//AttributeSet/Attribute"):
                            t_type = t.findtext("Type")
                            t_value = t.findtext("Value")
                            if t_type == "HGVS cDNA":
                                hgvs_c.append(t_value)
                            elif t_type == "HGVS protein":
                                hgvs_p.append(t_value)

                    # Filter if user entered a variant
                    if variant:
                        if not any(variant in s for s in hgvs_c):
                            continue  # skip this record if cDNA doesn't match

                    results.append({
                        "RCV ID": rcv_id or "N/A",
                        "Clinical Significance": significance or "N/A",
                        "Condition": condition or "N/A",
                        "Review Status": review_status or "N/A",
                        "HGVS cDNA": ", ".join(hgvs_c) if hgvs_c else "N/A",
                        "HGVS Protein": ", ".join(hgvs_p) if hgvs_p else "N/A"
                    })

                if results:
                    df = pd.DataFrame(results)
                    st.success(f"Found {len(df)} matching record(s).")
                    st.dataframe(df)
                else:
                    st.warning("No records match the HGVS variant entered.")

        except Exception as e:
            st.error(f"An error occurred: {e}")
