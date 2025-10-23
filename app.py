#app.py

import streamlit as st
import requests
import pandas as pd

st.set_page_config(page_title="ClinQuery", layout="centered")
st.title("ClinQuery")
st.markdown("Quick lookup: ClinVar + gnomAD (prototype)")

st.sidebar.header("Lookup")
query = st.sidebar.text_input("Enter gene name, HGVS c. notation, or rsID", "")

if st.sidebar.button("Search") and query.strip():
    st.info(f"Searching for: **{query}** — this is a prototype demo.")
    # Minimal placeholder logic: show the query back and a fake example table.
    # Replace this block later with real API calls to ClinVar/gnomAD.
    example = {
        "variant": [query],
        "clinvar_significance": ["Not queried yet (prototype)"],
        "gnomad_af": ["Not queried yet (prototype)"],
        "notes": ["Replace with real API results"]
    }
    df = pd.DataFrame(example)
    st.table(df)
else:
    st.write("Type a variant (e.g., `NPC1 c.3220A>T` or `rs12345`) in the sidebar and click **Search**.")
    st.markdown("When you're ready I can help add real ClinVar / gnomAD API calls.")
