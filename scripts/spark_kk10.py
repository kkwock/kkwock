import pandas as pd
import os
import math

# input name
SparkID = "014_kk10"

# Other Parameters
volume = 40 # Total volume in wells of Spark Plate
dilutionFactor = 5 # put in 5 for 5X dilution, etc.

# find file within local directories
filepath = "/Users/"

for root, dirs, files in os.walk(filepath):
    for file in files:
        if (SparkID + ".xlsx") in file:
            sparkFile = (os.path.join(root, file))

spark_data = pd.read_excel(sparkFile)

# function to pull out 260/320 abs and calculate concentrations
def colAbs(raw_spark, vol, dil):
    index260 = []
    index320 = []
    wells = []
    conc = []
    nmConc = []
    actualnM = []

    for names in raw_spark.iloc[35:(35 + 48), 0]:
        wells.append(names)
    for n260 in raw_spark.iloc[146:(146 + 48), 1]:
        index260.append(n260)
    for n320 in raw_spark.iloc[368:(368 + 48), 1]:
        index320.append(n320)
    # create one dataframe
    df = {"Wells": wells,
          "abs260": index260,
          "abs320": index320}
    df = pd.DataFrame(df)
    df = df.set_index("Wells")

    # Calculate concentration ng/uL ((260-320)*pi*4.69^2)/.002*4*vol)-8.5
    for i in range(len(df)):
        concalc = ((df.iloc[i, 0] - df.iloc[i, 1]) * math.pi * 4.69 ** 2 / (0.008 * vol)) - 8.5
        conc.append(round(concalc, 3))

    df.insert(2, "Concentration ng/uL", conc)

    # Calculate nM (Concentration/(660*65))*10^6
    for i in range(len(df)):
        nMole = ((df.iloc[i, 2]) / (660 * 65)) * 10 ** 6
        nmConc.append(round(nMole, 3))

    df.insert(3, "nM", nmConc)

    # Calculate "Actual nM in Plate"
    for i in range(len(df)):
        platenm = ((df.iloc[i, 3]) * dil)
        actualnM.append(round(platenm, 3))

    df.insert(4, "Actual nM in Plate", actualnM)

    # Average nM
    avenM = round(sum(df.iloc[:len(df), 4]) / len(df), 2)
    df.insert(5, "Average nM", avenM)

    print(df)

# write out a csv
colAbs(spark_data, volume, dilutionFactor)

