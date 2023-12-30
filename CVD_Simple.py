#!/usr/bin/env python3                                                          
import pandas as pd
import numpy as np
import sys

# Read the CSV file

import pandas as pd
import numpy as np
import sys

# Read the CSV file
data = pd.read_csv(sys.argv[1], encoding="ISO-8859-1")

# Define age and cholesterol points dictionaries
age_points = {
    (20, 34): -7, (35, 39): -3, (40, 44): 0, (45, 49): 3,
    (50, 54): 6, (55, 59): 8, (60, 64): 10, (65, 69): 12,
    (70, 74): 14, (75, 150): 16
}

chol_points = {
    (20, 39): {160: 0, 200: 4, 240: 8, 280: 11},
    (40, 49): {160: 0, 200: 3, 240: 6, 280: 8},
    (50, 59): {160: 0, 200: 2, 240: 4, 280: 5},
    (60, 69): {160: 0, 200: 1, 240: 2, 280: 3},
    (70, 79): {160: 0, 200: 1, 240: 1, 280: 2}
}

age_points_men = {
    (20, 34): -9, (35, 39): -4, (40, 44): 0, (45, 49): 3,
    (50, 54): 6, (55, 59): 8, (60, 64): 10, (65, 69): 11,
    (70, 74): 12, (75, 150): 13
}

chol_points_men = {
    (20, 39): {160: 0, 200: 4, 240: 7, 280: 11},
    (40, 49): {160: 0, 200: 3, 240: 5, 280: 8},
    (50, 59): {160: 0, 200: 2, 240: 3, 280: 5},
    (60, 69): {160: 0, 200: 1, 240: 1, 280: 3},
    (70, 79): {160: 0, 200: 0, 240: 0, 280: 1}
}

# Define a function to calculate CVD risk
def calculate_cvd_risk(row):
    age = row['Age']
    total_chol = row['Chol']
    hdl_chol = row['HDL']
    systolic_bp = row['BPsystolic']
    smoking = row['Cigarette']
    gender = row['Gender']

    # Determine which age and cholesterol points to use based on gender
    if gender == 'Female':
        age_points_dict = age_points
        chol_points_dict = chol_points
        bp_points_dict = {
            'Under 120': 0, '120-129': 1, '130-139': 2, '140-159': 3, '160 or higher': 4
        }
        smoking_points_dict = {
            'Age 20–39 years': 9, 'Age 40–49 years': 7, 'Age 50–59 years': 4,
            'Age 60–69 years': 2, 'Age 70–79 years': 1
        }
    else:  # Assume Male if not Female
        age_points_dict = age_points_men
        chol_points_dict = chol_points_men
        bp_points_dict = {
            'Under 120': 0, '120-129': 0, '130-139': 1, '140-159': 1, '160 or higher': 2
        }
        smoking_points_dict = {
            'Age 20–39 years': 8, 'Age 40–49 years': 5, 'Age 50–59 years': 3,
            'Age 60–69 years': 1, 'Age 70–79 years': 1
        }

    # Calculate age points
    age_group = next((group for group, points in age_points_dict.items() if group[0] <= age <= group[1]), (75, 150))
    age_points_value = age_points_dict.get(age_group)

    # Calculate total cholesterol points
    chol_group = chol_points_dict.get(age_group, {160: 0, 200: 1, 240: 1, 280: 2})
    chol_points_value = next((points for chol, points in chol_group.items() if total_chol >= chol), 2)

    # Calculate points for smoking status
    smoking_points = smoking_points_dict.get(age_group, 0) if smoking == 1 else 0

    # Calculate points for HDL cholesterol
    hdl_points = -1 if hdl_chol >= 60 else (0 if 50 <= hdl_chol < 60 else 1 if 40 <= hdl_chol < 50 else 2)

    # Calculate points for systolic blood pressure
    bp_group = next((group for group, points in bp_points_dict.items() if group == systolic_bp), 0)
    bp_points = bp_points_dict.get(bp_group, 0)

    # Calculate total points
    total_points = age_points_value + chol_points_value + smoking_points + hdl_points + bp_points

    # Calculate 10-year risk
    risk_percent = 1 if total_points < 0 else (total_points - (0 if gender == 'Male' else 9)) * 1 + 1

    return risk_percent

# Apply the function to the DataFrame
data['cvd_risk'] = data.apply(calculate_cvd_risk, axis=1)

# Save the results to a CSV file
data.to_csv('cvd_risk_results_male_female.csv', index=False)

print("Results for both genders saved to 'cvd_risk_results_male_female.csv'")

