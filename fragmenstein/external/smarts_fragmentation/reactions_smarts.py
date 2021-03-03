'''
Taken from https://github.com/xchem/opentrons/blob/master/xchem_ot/utils/reactions.py


FRANK's Knime

https://teams.microsoft.com/_#/school/files/General?threadId=19:eb7beb85d3c2403a90a83adf8b42453a@thread.tacv2&ctx=channel&context=Poised%2520design%2520-%2520Oakley%25202015&rootfolder=%252Fsites%252FCMD-XChemJointGroup%252FShared%2520Documents%252FGeneral%252FCode%252FPoised%2520design%2520-%2520Oakley%25202015

'''

REACTION_SMARTS = {

    'Amides': ['[#7:1][C;x0:2]=[O:3]>>[#7:1].Cl[C:2]=[O:3]'],
    'Benzimidazole': [
        '[#7:1]1[#6:9][#7:2][#6:4]2[#6:3]1[#6:8][#6:7][#6:6][#6:5]2>>[#7:1][c:3]1[c:8][c:7][c:6][c:5][c:4]1[#7:2].Cl[#6:9]=O'],
    'Benzoxazole': [
        '[#7:1]1[#6:9][#8:2][#6:4]2[#6:3]1[#6:8][#6:7][#6:6][#6:5]2>>[#7:1][c:3]1[c:8][c:7][c:6][c:5][c:4]1[#8:2].Cl[#6:9]=O'],
    'Ester_Coupling': ['[#8:1][C;x0:2]=[O:3]>>[#8:1].Cl[C:2]=[O:3]'],
    'Ether_Coupling': ['[CH0:2]-[#8R0:1]>>[#8:1].[#6:2]Br',
                       '[CH1:2]-[#8R0:1]>>[#8:1].[#6:2]Br',
                       '[CH2R0:2]-[#8:1]>>[#8:1].[#6:2]Br',
                       '[CH1R0:2]-[#8:1]>>[#8:1].[#6:2]Br',
                       '[CH2:2]-[#8R0:1]>>[#8:1].[#6:2]Br',
                       '[CH3:2]-[#8R0:1]>>[#8:1].[#6:2]Br',
                       '[CH0R0:2]-[#8:1]>>[#8:1].[#6:2]Br',
                       '[n:3][c:2]-[#8:1]>>[#8:1].[n:3][#6:2]Br'],
    'Indole': [
        '[c:10]1[c:9][nH:1][c:3]2[c:8]1[c:7][c:6][c:5][c:4]2>>N[N:1][c:3]1[c:8]([2H])[c:7][c:6][c:5][c:4]1.[C:10][C:9]=O'],
    'N-Alkylation': ['[CH0:2]-[#7R0:1]>>[#7:1].[#6:2]Br',
                     '[CH1:2]-[#7R0:1]>>[#7:1].[#6:2]Br',
                     '[CH2R0:2]-[#7:1]>>[#7:1].[#6:2]Br',
                     '[CH1R0:2]-[#7:1]>>[#7:1].[#6:2]Br',
                     '[CH2:2]-[#7R0:1]>>[#7:1].[#6:2]Br',
                     '[CH3:2]-[#7R0:1]>>[#7:1].[#6:2]Br',
                     '[CH0R0:2]-[#7:1]>>[#7:1].[#6:2]Br'],
    'Oxadiazole': ['[#6:6][c:4]1[n:5][o:3][c:1][n:2]1>>[O:2]-[C:1]=[O:3].[#6:6][C:4]#[N:5]'],
    'Reductive_Amination': ['[CH1:2]-[#7R0:1]>>[#7:1].[#6:2]=O',
                            '[CH2R0:2]-[#7:1]>>[#7:1].[#6:2]=O',
                            '[CH1R0:2]-[#7:1]>>[#7:1].[#6:2]=O',
                            '[CH2:2]-[#7R0:1]>>[#7:1].[#6:2]=O',
                            '[CH3:2]-[#7R0:1]>>[#7:1].[#6:2]=O'],
    'SNAr': ['[c:2][N:1][#6:3]>>[#6:3]-[#7:1].[c:2]Br'],
    'Sonogashira': ['[#6;a:1][C:2]#[C:3]>>[#6;a:1]Br.[C:2]#[C:3]'],
    'Sulfonamide': ['[#7:1][S:2](=[O:3])=[O:4]>>[#7:1].Cl[S:2](=[O:3])=[O:4]'],
    'Suzuki_Coupling': ['[#6;a:1]-[#6;a:2]>>[#6;a:2]Br.[#6;a:1]-[#5](-[#8])-[#8]'],
    'Triazole': ['[#6:6][n:1]1[c:3][c:4][n:2][n:5]1>>[#6:6][#7:1][#7:5]=[#7:2].[C:3]#[C:4]'],
    'Urea': ['[#7:1][C;x0]([#7:2])=O>>[#7:1].[#7:2]']

}


# REACTION_SMARTS = {
#     "amides": "[$([#1,*])]C(=O)N([$([#1,*])])[$([#1,*])]>>ClC([$([#1,*])])=O.[$([#1,*])]N[$([#1,*])]",
#     "sulfonamides": "[$([#1,*])]S(=O)(=O)N([$([#1,*])])[$([#1,*])]>>ClS([$([#1,*])])(=O)=O.[$([#1,*])]N[$([#1,*])]",
#     "sonogashira": "[$([#1,*])]C#C[$([#1,*])]>>[H]C#C[$([#1,*])].Br[$([#1,*])]",
#     "ureas": "[$([#1,*])]N([$([#1,*])])C(=O)N([$([#1,*])])[$([#1,*])]>>[$([#1,*])]N[$([#1,*])].[$([#1,*])]N[$([#1,*])]",
#     "arN": "[$([#1,*])]:c(:[$([#1,*])])N([$([#1,*])])[$([#1,*])]>>[Fe]N([$([#1,*])])[$([#1,*])].Brc(:[$([#1,*])]):[$([#1,*])]",
#     "suzuki": "[$([#1,*])]:c(:[$([#1,*])])-c(:[$([#1,*])]):[$([#1,*])]>>Brc(:[$([#1,*])]):[$([#1,*])].OB(O)c(:[$([#1,*])]):[$([#1,*])]"
#
# }

