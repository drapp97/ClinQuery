import streamlit as st
from Bio import Entrez
import pandas as pd
import xml.etree.ElementTree as ET

st.title("🧬 ClinQuery")
st.write("Search ClinVar for variant information quickly and easily.")

# Set your email for NCBI Entrez
Entrez.email = "daniel.rappaport@sickkids.ca"  # <-- replace with your email

query = st.text_input("Enter gene name or variant (e.g., NPC1 c.3044G>T):")

if query:
    st.info(f"Searching ClinVar for: **{query}** ...")
    try:
        # Step 1: Search ClinVar for matching records
        handle = Entrez.esearch(db="clinvar", term=query, retmax=5)
        record = Entrez.read(handle)
        handle.close()
        id_list = record["IdList"]

        if not id_list:
            st.warning("No results found.")
        else:
            # Step 2: Fetch detailed records
            handle = Entrez.efetch(db="clinvar", id=",".join(id_list), retmode="xml")
            data = handle.read()
            handle.close()

            root = ET.fromstring(data)
            results = []

            # Parse each ClinVar record
            for clinvar_record in root.findall(".//ClinVarSet"):
                variation_id = clinvar_record.findtext("ReferenceClinVarAssertion/ClinVarAccession/@Acc")
                significance = clinvar_record.findtext("ReferenceClinVarAssertion/ClinicalSignificance/Description")
                gene = clinvar_record.findtext("ReferenceClinVarAssertion/MeasureSet/Measure/MeasureRelationship/Target/Symbol")
                condition = clinvar_record.findtext("ReferenceClinVarAssertion/TraitSet/Trait/Name/ElementValue")
                review_status = clinvar_record.findtext("ReferenceClinVarAssertion/ClinicalSignificance/ReviewStatus")

                results.append({
                    "Variation ID": variation_id or "N/A",
                    "Clinical Significance": significance or "N/A",
                    "Gene": gene or "N/A",
                    "Condition": condition or "N/A",
                    "Review Status": review_status or "N/A"
                })

            df = pd.DataFrame(results)
            st.success(f"Found {len(df)} matching record(s).")
            st.dataframe(df)

    except Exception as e:
        st.error(f"An error occurred: {e}")
