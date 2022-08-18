import pandas as pd
import numpy as np
from os import listdir
import argparse

def find_nearest(array, value):
    #https://stackoverflow.com/questions/2566412/find-nearest-value-in-numpy-array
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx, array[idx]

def get_nearest_rows(df_lcms, times, indices, rt, mz):
    # find the nearest time from scans file to negID's rt
    selected_id, selected_rt = find_nearest(times, rt)
    # get the corresponding rows from the scans file for that time
    start = indices[selected_id]
    end = indices[selected_id+1] if selected_id + 1 < indices.size else len(df_lcms.index)
    return df_lcms.iloc[start:end, :].copy()

def filter_1(df, mz):
    # filtering could be applied but was generally not necessary
    selected_mz_ID = df['m/z'].sub(mz).abs().idxmin()
    Lowest_Intensity = df.at[selected_mz_ID,"Intensity"]*0.01
    df.drop(df[(df.Intensity < Lowest_Intensity) & (df.Zoom == 0)].index, inplace = True)
    # alternate approach
    # return df[ (df["Intensity"] > Lowest_Intensity) | (df["Zoom"] == 1) ]

def construct_table(df_lcms, lcms_times, lcms_indices, rt, fID, mz, dDa, filename):
    df = get_nearest_rows(df_lcms, lcms_times, lcms_indices, rt, mz)
    df.columns = ["Feature", "m/z", "Intensity"]
    df["Feature"] = fID
    df["Isotope"] = ""
    df["Formula"] = ""
    df["PredAbundance"] = ""
    df["File"] = filename
    df["Zoom"] = (np.abs(df["m/z"] - mz) < dDa) + 0
    filter_1(df, mz)
    return df.round({"m/z":5, "Intensity":1})

def extract_lcms_at_negID_rts(fn_output, df_NegID, fn_lcms, mz_id, rt_id, dDa):
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

        df = construct_table(df_lcms, lcms_times, lcms_indices, rt, fID, mz, dDa, filename)
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

def extract_all_lcms_scans(fn_negID, path, fn_output, mz_id, rt_id, dDa):
    '''
        load negID file
        get lcms csv files (corresponding to mzxml files)
        write to a target file for each scan
    '''
    df_NegID = pd.read_csv(fn_negID)
    fns_lcms = get_lcms_filenames(path)
    write_header(fn_output)
    for fn_lcms in fns_lcms:
        extract_lcms_at_negID_rts(fn_output, df_NegID, fn_lcms, mz_id, rt_id, dDa)

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

def run_program(*args, **kwargs):
    print(locals())
    extract_all_lcms_scans(fn_negID = kwargs["fn_negID"], path = kwargs["path"], fn_output = kwargs["fn_output"],
         mz_id = kwargs["mz_id"], rt_id = kwargs["rt_id"], dDa = kwargs["dDa"])

if __name__ == '__main__':
    args = parser.parse_args()
    run_program(**vars(args))