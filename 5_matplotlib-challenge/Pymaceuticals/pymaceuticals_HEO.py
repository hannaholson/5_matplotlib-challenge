# Dependencies and Setup
import matplotlib.pyplot as plt
import pandas as pd
import scipy.stats as st

# Study data files
mouse_metadata_path = "data/Mouse_metadata.csv"
study_results_path = "data/Study_results.csv"

# Read the mouse data and the study results
mouse_metadata = pd.read_csv(mouse_metadata_path)
study_results = pd.read_csv(study_results_path)

# Combine the data into a single dataset
# Merge LEFT
full_study_data = pd.merge(study_results, mouse_metadata, how = "left", on = "Mouse ID")

# Display the data table for preview
full_study_data

# Checking the number of mice.
# Using .unique

len(full_study_data["Mouse ID"].unique())

# Getting the duplicate mice by ID number that shows up for Mouse ID and Timepoint. 
# Using .loc and .unique

mouse_id_dupes = full_study_data.loc[full_study_data.duplicated(subset=['Mouse ID', 'Timepoint']), 'Mouse ID'].unique()

# Display the data table for preview

mouse_id_dupes

# Optional: Get all the data for the duplicate mouse ID. 
# Using .loc

mouse_data_dupes = full_study_data.loc[full_study_data['Mouse ID'] == 'g989']

# Display the data table for preview

mouse_data_dupes

# Create a clean DataFrame by dropping the duplicate mouse by its ID.
# Using .isin

full_study_data_clean = full_study_data[full_study_data['Mouse ID'].isin(mouse_id_dupes) == False]

# Display the data table for preview

full_study_data_clean

# Checking the number of mice in the clean DataFrame.
# Using .unique

len(full_study_data_clean['Mouse ID'].unique())

# Generate a summary statistics table of mean, median, variance, standard deviation, and SEM of the tumor volume for each regimen

# Use groupby and summary statistical methods to calculate the following properties of each drug regimen:
    # mean, median, variance, standard deviation, and SEM of the tumor volume. 

means = full_study_data_clean.groupby('Drug Regimen').mean()['Tumor Volume (mm3)']
medians = full_study_data_clean.groupby('Drug Regimen').median()['Tumor Volume (mm3)']
variances = full_study_data_clean.groupby('Drug Regimen').var()['Tumor Volume (mm3)']
sds = full_study_data_clean.groupby('Drug Regimen').std()['Tumor Volume (mm3)']
sems = full_study_data_clean.groupby('Drug Regimen').sem()['Tumor Volume (mm3)']

# Assemble the resulting series into a single summary dataframe.

summary_table = pd.DataFrame({"Mean Tumor Volume":means,
                              "Median Tumor Volume":medians,
                              "Tumor Volume Variance":variances,
                              "Tumor Volume Std. Dev.":sds,
                              "Tumor Volume Std. Err.":sems})

# Display the data table for preview

summary_table

# Generate a summary statistics table of mean, median, variance, standard deviation, and SEM of the tumor volume for each regimen

# Using the aggregation method, produce the same summary statistics in a single line

summary_table = full_study_data_clean.groupby("Drug Regimen").agg({"Tumor Volume (mm3)":["mean","median","var","std","sem"]})

# Display the data table for preview

summary_table

# Generate a bar plot showing the total number of measurements taken on each drug regimen using pandas.
# Using .plot(kind="bar")

counts = full_study_data_clean['Drug Regimen'].value_counts()
counts.plot(kind="bar")
plt.xlabel("Drug Regimen")
plt.xticks(rotation=90)
plt.ylabel("Total Number of Measurements")

# Display the bar graph for preview
plt.show()

# Generate a bar plot showing the total number of measurements taken on each drug regimen using pyplot.
# Using .plot(kind="bar")

counts = full_study_data_clean['Drug Regimen'].value_counts()
plt.bar(counts.index.values,counts.values)
plt.xlabel("Drug Regimen")
plt.xticks(rotation=90)
plt.ylabel("Total Number of Measurements")


# Display the bar graph for preview
plt.show()

# Generate a pie plot showing the distribution of female versus male mice using pandas
# Using .plot(kind="pie")
counts = full_study_data_clean.Sex.value_counts()
counts.plot(kind="pie",autopct='%1.1f%%')

# Display the pie plot for preview
plt.show()

# Generate a pie plot showing the distribution of female versus male mice using pyplot
# Using .plot(kind="pie")

counts = full_study_data_clean.Sex.value_counts()
plt.pie(counts.values,labels=counts.index.values,autopct='%1.1f%%')
plt.ylabel('Sex')

# Display the pie plot for preview

plt.show()

# Calculate the final tumor volume of each mouse across four of the treatment regimens:  
# Capomulin, Ramicane, Infubinol, and Ceftamin
# Using .groupby and .max

# Start by getting the last (greatest) timepoint for each mouse

max_tumor_volume = full_study_data_clean.groupby(['Mouse ID'])['Timepoint'].max()
max_tumor_volume = max_tumor_volume.reset_index()

# Merge this group df with the original dataframe to get the tumor volume at the last timepoint
# Using .merge

merged_data = max_tumor_volume.merge(full_study_data_clean, on = ['Mouse ID', 'Timepoint'], how = "left")

# Put treatments into a list for for loop (and later for plot labels)

list_o_treatments = ['Capomulin', 'Ramicane', 'Infubinol', 'Ceftamin']

# Create empty list to fill with tumor vol data (for plotting)

list_o_tumor_volumes = []

# Calculate the IQR and quantitatively determine if there are any potential outliers. 

for drug in list_o_treatments:
    
    # Locate the rows which contain mice on each drug and get the tumor volumes
    # Using .loc
    
    tumor_volume_final = merged_data.loc[merged_data['Drug Regimen'] == drug, 'Tumor Volume (mm3)']
    
    # add subset 
    # Using .append
    
    list_o_tumor_volumes.append(tumor_volume_final)
    
    # Determine outliers using upper and lower bounds
    # Using .quantile
    
    quartiles = tumor_volume_final.quantile([.25,.5,.75])
    
    lowerq = quartiles[0.25]
    upperq = quartiles[0.75]
    iqr = upperq-lowerq
   
    lower_bound = lowerq - (1.5*iqr)
    upper_bound = upperq + (1.5*iqr)
    
    outliers = tumor_volume_final.loc[(tumor_volume_final < lower_bound) | (tumor_volume_final > upper_bound)]
    
    # Display Outliers
    print(f"{drug}'s outliers: {outliers}")

# Generate a box plot of the final tumor volume of each mouse across four regimens of interest
# Using .boxplot

blue_out = dict(markerfacecolor = 'green', markersize = 10)
plt.boxplot(list_o_tumor_volumes, labels = list_o_treatments, flierprops = blue_out)
plt.ylabel('Final Tumor Volume (mm3)')

# Display box plot

plt.show()

# Generate a line plot of tumor volume vs. time point for a mouse treated with Capomulin
# Using .plot

capomulin_table = full_study_data_clean.loc[full_study_data_clean['Drug Regimen'] == 'Capomulin']
mousedata = capomulin_table.loc[capomulin_table['Mouse ID'] == '1509']
plt.plot(mousedata['Timepoint'], mousedata['Tumor Volume (mm3)'])
plt.xlabel('Timepoint (days)')
plt.ylabel('Tumor Volume (mm3)')
plt.title('Capomulin treatment of mouse 1509')

# Display scatter plot

plt.show()

# Generate a scatter plot of average tumor volume vs. mouse weight for the Capomulin regimen
# Using .scatter

capomulin_table = full_study_data_clean.loc[full_study_data_clean['Drug Regimen'] == 'Capomulin']
capomulin_average = capomulin_table.groupby(['Mouse ID']).mean()
plt.scatter(capomulin_average['Weight (g)'], capomulin_average['Tumor Volume (mm3)'])
plt.xlabel('Weight (g)')
plt.ylabel('Average Tumor Volume (mm3)')

# Display scatter plot
plt.show

# Calculate the correlation coefficient and linear regression model 
# for mouse weight and average tumor volume for the Capomulin regimen

corr=round(st.pearsonr(capomulin_average['Weight (g)'],capomulin_average['Tumor Volume (mm3)'])[0],2)
print(f"The correlation between mouse weight and the average tumor volume is {corr}")
model = st.linregress(capomulin_average['Weight (g)'],capomulin_average['Tumor Volume (mm3)'])

y_values = capomulin_average['Weight (g)']*model[0]+model[1]
plt.scatter(capomulin_average['Weight (g)'],capomulin_average['Tumor Volume (mm3)'])
plt.plot(capomulin_average['Weight (g)'],y_values,color="red")
plt.xlabel('Weight (g)')
plt.ylabel('Average Tumor Volume (mm3)')

# Display scatter plot

plt.show()

# THE END