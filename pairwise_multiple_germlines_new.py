# Finds and scores all framework mutations from input antibody file (csv format). Outputs normalized FR scores. 
# Verbose mode prints full pairwise alignment of each antibody.
# output_mutations option creates a csv with all antibody scores

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from Bio import pairwise2
from Bio.pairwise2 import format_alignment

def get_position(pos, germ):
    """ Place gaps for IMGT numbering scheme """
    try:
        return positions[positions[germ] == pos].index[0]
    except:
        return 0

dict = {
    "1-2": "QVQLVQSGAEVKKPGASVKVSCKASGYTFTGYYMHWVRQAPGQGLEWMGRINPNSGGTNYAQKFQGRVTSTRDTSISTAYMELSRLRSDDTVVYYCAR",
    "1-3": "QVQLVQSGAEVKKPGASVKVSCKASGYTFTSYAMHWVRQAPGQRLEWMGWINAGNGNTKYSQKFQGRVTITRDTSASTAYMELSSLRSEDTAVYYCAR",
    "1-24": "QVQLVQSGAEVKKPGASVKVSCKVSGYTLTELSMHWVRQAPGKGLEWMGGFDPEDGETIYAQKFQGRVTMTEDTSTDTAYMELSSLRSEDTAVYYCAT",
    "1-46": "QVQLVQSGAEVKKPGASVKVSCKASGYTFTSYYMHWVRQAPGQGLEWMGIINPSGGSTSYAQKFQGRVTMTRDTSTSTVYMELSSLRSEDTAVYYCAR",
    "1-69": "QVQLVQSGAEVKKPGSSVKVSCKASGGTFSSYAISWVRQAPGQGLEWMGGIIPIFGTANYAQKFQGRVTITADESTSTAYMELSSLRSEDTAVYYCAR",
    "2-5": "QITLKESGPTLVKPTQTLTLTCTFSGFSLSTSGVGVGWIRQPPGKALEWLALIYWNDDKRYSPSLKSRLTITKDTSKNQVVLTMTNMDPVDTATYYCAHR",
    "3-7": "EVQLVESGGGLVQPGGSLRLSCAASGFTFSSYWMSWVRQAPGKGLEWVANIKQDGSEKYYVDSVKGRFTISRDNAKNSLYLQMNSLRAEDTAVYYCAR",
    "3-9": "EVQLVESGGGLVQPGRSLRLSCAASGFTFDDYAMHWVRQAPGKGLEWVSGISWNSGSIGYADSVKGRFTISRDNAKNSLYLQMNSLRAEDTALYYCAKD",
    "3-20": "EVQLVESGGGVVRPGGSLRLSCAASGFTFDDYGMSWVRQAPGKGLEWVSGINWNGGSTGYADSVKGRFTISRDNAKNSLYLQMNSLRAEDTALYHCAR",
    "3-21": "EVQLVESGGGLVKPGGSLRLSCAASGFTFSSYSMNWVRQAPGKGLEWVSSISSSSSYIYYADSVKGRFTISRDNAKNSLYLQMNSLRAEDTAVYYCAR",
    "3-23": "EVQLLESGGGLVQPGGSLRLSCAASGFTFSSYAMSWVRQAPGKGLEWVSAISGSGGSTYYADSVKGRFTISRDNSKNTLYLQMNSLRAEDTAVYYCAK",
    "3-30": "QVQLVESGGGVVQPGRSLRLSCAASGFTFSSYAMHWVRQAPGKGLEWVAVISYDGSNKYYADSVKGRFTISRDNSKNTLYLQMNSLRAEDTAVYYCAR",
    "3-33": "QVQLVESGGGVVQPGRSLRLSCAASGFTFSSYGMHWVRQAPGKGLEWVAVIWYDGSNKYYADSVKGRFTISRDNSKNTLYLQMNSLRAEDTAVYYCAR",
    "3-48": "EVQLVESGGGLVQPGGSLRLSCAASGFTFSSYSMNWVRQAPGKGLEWVSYISSSSSTIYYADSVKGRFTISRDNAKNSLYLQMNSLRAEDTAVYYCAR",
    "3-66": "EVQLVESGGGLVQPGGSLRLSCAASGFTVSSNYMSWVRQAPGKGLEWVSVIYSGGSTYYADSVKGRFTISRDNSKNTLYLQMNSLRAEDTAVYYCAR",
    "3-74": "EVQLVESGGGLVQPGGSLRLSCAASGFTFSSYWMHWVRQAPGKGLVWVSRINSDGSSTSYADSVKGRFTISRDNAKNTLYLQMNSLRAEDTAVYYCAR",
    "4-4": "QVQLQESGPGLVKPPGTLSLTCAVSGGSISSSNWWSWVRQPPGKGLEWIGEIYHSGSTNYNPSLKSRVTISVDKSKNQFSLKLSSVTAADTAVYCCAR",
    "4-30-4": "QVQLQESGPGLVKPSQTLSLTCTVSGGSISSGDYYWSWIRQPPGKGLEWIGYIYYSGSTYYNPSLKSRVTISVDTSKNQFSLKLSSVTAADTAVYYCAR",
    "4-31": "QVQLQESGPGLVKPSQTLSLTCTVSGGSISSGGYYWSWIRQHPGKGLEWIGYIYYSGSTYYNPSLKSLVTISVDTSKNQFSLKLSSVTAADTAVYYCAR",
    "4-34": "QVQLQQWGAGLLKPSETLSLTCAVYGGSFSGYYWSWIRQPPGKGLEWIGEINHSGSTNYNPSLKSRVTISVDTSKNQFSLKLSSVTAADTAVYYCAR",
    "4-39": "QLQLQESGPGLVKPSETLSLTCTVSGGSISSSSYYWGWIRQPPGKGLEWIGSIYYSGSTYYNPSLKSRVTISVDTSKNQFSLKLSSVTAADTAVYYCAR",
    "4-59": "QVQLQESGPGLVKPSETLSLTCTVSGGSISSYYWSWIRQPPGKGLEWIGYIYYSGSTNYNPSLKSRVTISVDTSKNQFSLKLSSVTAADTAVYYCAR",
    "4-61": "QVQLQESGPGLVKPSETLSLTCTVSGGSVSSGSYYWSWIRQPPGKGLEWIGYIYYSGSTNYNPSLKSRVTISVDTSKNQFSLKLSSVTAADTAVYYCAR",
    "5-51": "EVQLVQSGAEVKKPGESLKISCKGSGYSFTSYWIGWVRQMPGKGLEWMGIIYPGDSDTRYSPSFQGQVTISADKSISTAYLQWSSLKASDTAMYYCAR",
    "6-1": "QVQLQQSGPGLVKPSQTLSLTCAISGDSVSSNSAAWNWIRQSPSRGLEWLGRTYYRSKWYNDYAVSVKSRITINPDTSKNQFSLQLNSVTPEDTAVYYCAR"
}

verbose = False                              #True: output alignment and indinvidual scores, False: plot results
output_csv = True                            #output csv file of calculated scores
output_mutations =  True                     #output csv file with all antibody scores

#input files
ab_filename = "FDA_Abs.csv"                #input antibody file (use "Flagged" in name for phase identification)
norm_filename = "normalization.csv"        #input normalization constants
numbering = "IMGT_num.csv"                 #index to IMGT numbering scheme

#read input files
Abs = pd.read_csv(ab_filename)
norm = pd.read_csv(norm_filename, index_col=0)
num_sequences = len(Abs["Name"])
positions = pd.read_csv(numbering, index_col=0)

names, x, y, mut, scores, gm = [], [], [], [], [], []

#pairwise alignment with germline, calculate F scores
for seq in range(num_sequences):
    name = Abs["Name"][seq]
    Ab = Abs["VH"][seq]
    germline_gene = Abs["Germline VH"][seq].replace("VH", "")
    if verbose: print(name, germline_gene)
    filename = "pssm_" + germline_gene + ".csv"         #match germline to pssm
    germline = dict[germline_gene]
    gm_count = len(germline)
    data = pd.read_csv(filename, index_col=0).replace("#WT", np.nan).apply(pd.to_numeric)       #read pssm
    cols = list(data.columns)
    alignments = pairwise2.align.globalms(germline, Ab, 2, -1, -5, -0.1)  #generate paiwise alignment of germline/Ab
    aln = format_alignment(*alignments[0])
    if verbose: print(aln)
    start1 = [aln.find(" "), aln.find("|"), aln.find(".")]
    start = min(start1)         #find starting search position
    mutations = []
    ins_count = 0
    for i in range(0,start):
        if aln[i] == "-":
            ins_count += 1
        if aln[start + i] == ".":
            g = aln[i]
            m = aln[start * 2 + i]
            pos = get_position(i - ins_count + 1, "VH"+germline_gene)
            if str(pos) in cols and m != "-" and ~np.isnan(data.at[m,str(pos)]):
                if verbose: print(g, pos, m, " ", round(data.at[m,str(pos)],1), sep='')
                mut.append(str(g+str(pos)+m))                                             
                scores.append(round(data.at[m,str(pos)],3))                                 #read mutation score
                gm.append(germline_gene)                                                  
                mutations.append(data.at[m,str(pos)])                                       
    if verbose: print("Sum scores: ", round(sum(mutations),1))
    if verbose: print("Number Mutations: ", len(mutations))
    if verbose: print("\n")
    names.append(name)
    x.append(len(mutations))
    if len(mutations) != 0:
        y.append(sum(mutations)/(norm.at[germline_gene, "c0"]*len(mutations) + norm.at[germline_gene, "c1"]*len(mutations)**2))         #calculate FR score
    else:
        y.append(1)         #replace undefined scores with score of 1

#log in dataframe, output to csv
df = pd.DataFrame({"Num_mut": x, "FR": y, "Name": names})
df_mut = pd.DataFrame({"Mutation": mut, "Score": scores, "Germline": gm})
outputname = "FRscores_normalized_" + ab_filename
if output_csv: df.to_csv(outputname, index=False)
if output_mutations: df_mut.to_csv("new_mutations.csv", index=False)

#generate scatter plot with labels
if not verbose:
    plt.figure(figsize=(12,8))
    sns.set(font_scale=1.5)
    sns.set_style("whitegrid")
    sns.scatterplot(x="Num_mut", y="FR", data=df)
    plt.title(ab_filename.replace(".csv", "").replace("_", " "))
    plt.xlabel("Number of mutations")
    plt.ylim(0,1.6)
    plt.ylabel("FR Score")
    sns.set(font_scale=1)
    for i, txt in enumerate(names):
        plt.annotate(txt, (x[i], y[i]), textcoords="offset points", xytext=(0.2,2), ha="left")
    plt.savefig("FDA_Abs.svg")
    plt.show()
