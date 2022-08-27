import pandas as pd
import numpy as np
from os import listdir
import argparse

def get_df_iso():
    arr_iso_rmass = [0.99939,1.9958,1.00336,2.00671,1.99705,3.9941,5.99115,7.9882,2.00424,1.99795,3.99591,0.99704]
    arr_iso_symbol = ["33S","34S","13C","13C2","37Cl","37Cl2","37Cl3","37Cl4","18O","81Br","81Br2","15N"]
    arr_iso_rint = [0.7893,4.4306,1.0816,0.0117,32.3977,10.4961,3.4005,0.8501,0.2005,97.5114,48.7557,0.3613]

    df = pd.DataFrame({
        'symbol': arr_iso_symbol
        ,'rmass': arr_iso_rmass
        ,'intensity': arr_iso_rint
        })
    return df

def find_nearest(array, value):
    #https://stackoverflow.com/questions/2566412/find-nearest-value-in-numpy-array
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx, array[idx]

def get_nearest_rows(df_lcms, times, indices, rt):
    # find the nearest time from scans file to negID's rt
    selected_id, selected_rt = find_nearest(times, rt)
    # get the corresponding rows from the scans file for that time
    start = indices[selected_id]
    end = indices[selected_id+1] if selected_id + 1 < indices.size else len(df_lcms.index)
    return df_lcms.iloc[start:end, :].copy()

def remove_under_lowest_percentage(df, dmz, lowest_percentage):
    # filtering could be applied but was generally not necessary
    selected_mz_ID = dmz.idxmin()
    Lowest_Intensity = df.at[selected_mz_ID,"Intensity"]*lowest_percentage
    df.drop(df[(df.Intensity < Lowest_Intensity) & (df.Zoom == 0)].index, inplace = True)
    # alternate approach
    # return df[ (df["Intensity"] > Lowest_Intensity) | (df["Zoom"] == 1) ]

def iso_string_gen(df_iso_matches, match_i):
    ''' For a particular index, compose the isotopic string. '''
    iso = df_iso_matches.iloc[match_i]
    return f"{iso.symbol}({iso.intensity:0.2f}%;{iso.d_rmass_to_dmzscan_i:0.5f}Da)"

def process_isotopes(df, dmz, isotope_threshold):
    selected_mz_ID = dmz.idxmin()
    df.at[selected_mz_ID, "Isotope"] = "M"
    tmz_scan = df.at[selected_mz_ID,"m/z"]
    # gets the index of the last row in the zoom window
    last_zoom_idx = df[df["Zoom"] == 1].index[-1]
    dmz_scan = df.loc[selected_mz_ID+1:last_zoom_idx, "m/z"] - tmz_scan
    df_iso = get_df_iso()

    for dmz_i in dmz_scan.index.values:
        # lists isotopic matches where the relative mass is close enough to the distance from a m/z to the tmz_scan
        df_iso["d_rmass_to_dmzscan_i"] = np.abs(df_iso.rmass - dmz_scan[dmz_i])
        df_iso_matches = df_iso[ df_iso.d_rmass_to_dmzscan_i < isotope_threshold ].sort_values(by=["d_rmass_to_dmzscan_i", "intensity"])
        iso_strings = []
        for match_i in range(df_iso_matches.shape[0]):
            iso_strings += [iso_string_gen(df_iso_matches, match_i)]
        iso_string = ";".join(iso_strings)
        df.at[dmz_i, "Isotope"] = iso_string

def construct_table(df_lcms, lcms_times, lcms_indices, rt, fID, mz, dDa, filename, lowest_percetange, isotope_threshold):
    df = get_nearest_rows(df_lcms, lcms_times, lcms_indices, rt)
    df.columns = ["Feature", "m/z", "Intensity"]
    df["Feature"] = fID
    df["Isotope"] = ""
    df["Formula"] = ""
    df["PredAbundance"] = ""
    df["File"] = filename
    dmz = np.abs(df["m/z"] - mz)
    df["Zoom"] = (dmz < dDa) + 0
    if lowest_percetange > 0:
        remove_under_lowest_percentage(df, dmz, lowest_percetange)
    if isotope_threshold > 0:
        process_isotopes(df, dmz, isotope_threshold)
    return df.round({"m/z":5, "Intensity":1})

def extract_lcms_at_negID_rts(fn_output, df_NegID, fn_lcms, mz_id, rt_id, dDa, lowest_percentage, isotope_threshold):
    '''
        load lcms scan
        get unique times
        for each row in negID
        get rows nearest to that retention time from lcms scan
        print the rows to the csv
    '''
    df_lcms = pd.read_csv(fn_lcms)
    lcms_times, lcms_indices = np.unique(df_lcms.rt, return_index=True)
    filename = fn_lcms.split("/")[-1]

    for row_id in range(len(df_NegID.index)):
        mz = df_NegID.iat[row_id,mz_id]      #mass to charge ratio
        rt = df_NegID.iat[row_id,rt_id] * 60 #retention time
        fID = df_NegID.at[row_id,"row.ID"]   #feature number

        df = construct_table(df_lcms, lcms_times, lcms_indices, rt, fID, mz, dDa, filename, lowest_percentage, isotope_threshold)
        df.to_csv(fn_output, mode='a', index=False, header=False)

            # example:
            # 227492,151.07624,33665.4,,,,FoamSample01_Targeted_Neg.csv,0
            # 227492,152.07179,39222.5,,,,FoamSample01_Targeted_Neg.csv,0

def get_lcms_filenames(path):
    '''
        get lcms csv files (corresponding to mzxml files)
        filenames for liquid chromatography mass spectrometry (LCMS) scans
    '''
    # load csv scan files
    allfiles = listdir(path)
    mzxmlfiles = sorted([i for i in allfiles if i[-5:].lower() == "mzxml"])
    # csvfiles = [i.replace("mzXML", "csv") for i in mzxmlfiles]
    csvfiles = [i[:-5] + "csv" for i in mzxmlfiles] # catches both mzXML and mzxml
    # "FoamSample01_Targeted_Neg.mzxml"[:-5] == "FoamSample01_Targeted_Neg."
    # +"csv" = "FoamSample01_Targeted_Neg.csv"
    csvfiles = [path+i for i in csvfiles if i in allfiles]
    return csvfiles

def write_header(fn_output):
    cols = ["Feature", "m/z", "Intensity", "Isotope", "Formula", "PredAbundance", "File", "Zoom"]
    df = pd.DataFrame({i:[] for i in cols})
    df.to_csv(fn_output, index=False)

def extract_all_lcms_scans(fn_negID, path, fn_output, mz_id, rt_id, dDa, lowest_percentage, isotope_threshold):
    '''
        load negID file
        get lcms csv files (corresponding to mzxml files)
        write to a target file for each scan
    '''
    df_NegID = pd.read_csv(fn_negID)
    fns_lcms = get_lcms_filenames(path)
    write_header(fn_output)
    for fn_lcms in fns_lcms:
        extract_lcms_at_negID_rts(fn_output, df_NegID, fn_lcms, mz_id, rt_id, dDa, lowest_percentage, isotope_threshold)

parser = argparse.ArgumentParser(description='msspectra')

parser.add_argument('--fn_negID',
                    dest='fn_negID',
                    type=str,
                    help='NegIDed file path',
                    default="NegIDed_FIN.csv")
parser.add_argument('--path',
                    dest='path',
                    type=str,
                    help='Test files path',
                    default="./")
parser.add_argument('--fn_output',
                    dest='fn_output',
                    type=str,
                    help='output file path',
                    default="scan_extractions_output.csv")
parser.add_argument('--mz_id',
                    dest='mz_id',
                    type=int,
                    help='NegID column index (Mass to charge ratio)',
                    default=5)
parser.add_argument('--rt_id',
                    dest='rt_id',
                    type=int,
                    help='NegID column index (Retention time)',
                    default=6)
parser.add_argument('--dDa',
                    dest='dDa',
                    type=int,
                    help='Zoom window in Daltons',
                    default=5)
parser.add_argument('--filter',
                    dest='lowest_percentage',
                    type=float,
                    help='Filtering cutoff value, ex: 0.01 is 1 percent of target intensity.',
                    default=0)
parser.add_argument('--isotope_threshold',
                    dest='isotope_threshold',
                    type=float,
                    help='Isotope cutoff value for determining a match, ex: 0.005Da.',
                    default=0.005)

def run_program(*args, **kwargs):
    print(locals())
    extract_all_lcms_scans(fn_negID = kwargs["fn_negID"], path = kwargs["path"], fn_output = kwargs["fn_output"],
         mz_id = kwargs["mz_id"], rt_id = kwargs["rt_id"], dDa = kwargs["dDa"]
         , lowest_percentage = kwargs["lowest_percentage"], isotope_threshold = kwargs["isotope_threshold"])

if __name__ == '__main__':
    args = parser.parse_args()
    run_program(**vars(args))