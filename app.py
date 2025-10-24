# app.py
import streamlit as st
import requests
import pandas as pd
import urllib.parse

st.set_page_config(page_title="ClinQuery", layout="wide")
st.title("🧬 ClinQuery")
st.write("Search ClinVar variants (via NLM ClinicalTables API). Enter a gene and optionally an HGVS cDNA substring to filter results.")

# Sidebar inputs
st.sidebar.header("Variant lookup")
gene = st.sidebar.text_input("Gene symbol (e.g., BRCA1)", "")
hgvs_filter = st.sidebar.text_input("Optional: HGVS cDNA substring (e.g., c.68_69delAG)", "")
max_results = st.sidebar.number_input("Max results to request (max 500)", min_value=1, max_value=500, value=50, step=1)

if st.sidebar.button("Search"):
    gene = gene.strip()
    if not gene:
        st.sidebar.error("Please enter a gene symbol.")
    else:
        st.info(f"Searching for gene: **{gene}** (assembly GRCh37/na dataset) ...")

        # Build ClinicalTables API request
        base = "https://clinicaltables.nlm.nih.gov/api/variants/v4/search"
        params = {
            "terms": gene,
            "count": max_results,   # page size
            "offset": 0,
            # df = display fields we want in the display output (human readable)
            "df": "VariationID,HGVS_c,HGVS_p,GeneSymbol,dbSNP,phenotype,GenomicLocation",
            # ef = extra fields we want in the extra-data hash (json)
            "ef": "GeneSymbol,HGVS_c,HGVS_p,VariationID,dbSNP,phenotype,GenomicLocation"
        }

        url = base + "?" + urllib.parse.urlencode(params, safe=',:')
        try:
            r = requests.get(url, timeout=20)
            r.raise_for_status()
            resp = r.json()
            # Response format (see docs):
            # resp[0] = total number available (int)
            # resp[1] = list of codes (VariationIDs)
            # resp[2] = dict of extra fields (each key -> list of values aligned with codes)
            # resp[3] = display strings (array per code)
            total = resp[0] if len(resp) > 0 else 0
            codes = resp[1] if len(resp) > 1 else []
            extra = resp[2] if len(resp) > 2 else {}
            display = resp[3] if len(resp) > 3 else []

            if not codes:
                st.warning(f"No variants found for gene {gene}. (Try increasing Max results.)")
            else:
                # Build dataframe from 'codes' + extra fields
                df_dict = {
                    "VariationID": codes
                }
                # include selected extra fields if present
                for key in ["GeneSymbol", "HGVS_c", "HGVS_p", "dbSNP", "phenotype", "GenomicLocation"]:
                    df_dict[key] = extra.get(key, ["N/A"] * len(codes))

                df = pd.DataFrame(df_dict)

                # Normalize columns
                if "HGVS_c" in df.columns:
                    df["HGVS_c"] = df["HGVS_c"].astype(str).replace("None", "")
                else:
                    df["HGVS_c"] = ""

                # Optional filter by user-entered HGVS substring (case-insensitive)
                if hgvs_filter.strip():
                    pat = hgvs_filter.strip()
                    mask = df["HGVS_c"].str.contains(pat, case=False, na=False)
                    df_filtered = df[mask]
                    if df_filtered.empty:
                        st.warning(f"No records matched HGVS substring '{hgvs_filter}'. Showing unfiltered results for gene {gene}.")
                        display_df = df
                    else:
                        st.success(f"Found {len(df_filtered)} record(s) matching HGVS substring '{hgvs_filter}'.")
                        display_df = df_filtered
                else:
                    st.success(f"Found {len(df)} variant(s) (showing up to {max_results}).")
                    display_df = df

                # Add ClinVar clickable link column
                def clinvar_link(vid):
                    try:
                        return f"https://www.ncbi.nlm.nih.gov/clinvar/variation/{vid}/"
                    except:
                        return ""

                display_df = display_df.reset_index(drop=True)
                display_df["ClinVar Link"] = display_df["VariationID"].apply(lambda v: clinvar_link(v))

                # Show a table with clickable links
                # Streamlit can't make a DataFrame column directly clickable; use st.write with markdown in one column
                display_df_for_show = display_df.copy()
                display_df_for_show["ClinVar Link"] = display_df_for_show["ClinVar Link"].apply(lambda u: f"[link]({u})" if u else "")
                st.dataframe(display_df_for_show, use_container_width=True)

                # Also provide CSV download
                csv = display_df.to_csv(index=False)
                st.download_button("Download results as CSV", csv, file_name=f"{gene}_clinvar_variants.csv", mime="text/csv")

        except requests.exceptions.RequestException as e:
            st.error(f"Network/API error: {e}")
        except Exception as e:
            st.error(f"Unexpected error: {e}")
