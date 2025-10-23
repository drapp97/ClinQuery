import streamlit as st
import requests
import pandas as pd

st.title("🧬 ClinQuery")
st.write("Search ClinVar for variant information quickly and easily.")

# Search input
query = st.text_input("Enter variant or gene name (e.g., BRCA1 c.68_69delAG):")

if query:
    st.info(f"Searching ClinVar for: **{query}** ...")
    try:
        # Query the ClinVar API
        url = f"https://api.ncbi.nlm.nih.gov/variation/v0/clinvar?term={query}"
        response = requests.get(url)
        response.raise_for_status()
        data = response.json()

        # Extract relevant information
        if "data" in data and data["data"]:
            results = []
            for item in data["data"]:
                record = {
                    "Variation ID": item.get("variation_id", "N/A"),
                    "Clinical Significance": item.get("clinical_significance", {}).get("description", "N/A"),
                    "Gene": item.get("gene", {}).get("symbol", "N/A"),
                    "Condition": item.get("condition", [{}])[0].get("name", "N/A"),
                    "Review Status": item.get("review_status", "N/A"),
                    "Last Updated": item.get("last_updated", "N/A")
                }
                results.append(record)
            
            df = pd.DataFrame(results)
            st.success(f"Found {len(df)} matching record(s).")
            st.dataframe(df)
        else:
            st.warning("No results found for that query.")
    except Exception as e:
        st.error(f"An error occurred: {e}")
