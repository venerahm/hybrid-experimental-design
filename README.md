# repository-hybrid-trials-QKWHF6125W

### sim_data_for_proximal.R

This file contains simulated data for the proximal analysis. It includes responders and nonresponders. The variables included are detailed below:

- **id**: participant id
- **female**: binary sex (1: female; 0: male) 
- **sex_centered**: binary sex (as defined above) mean centered 
- **bmi_bl**: body mass index (bmi) at baseline in kg/m^2
- **bmi_bl_centered**: bmi mean centered
- **Z1**: initial intervention (1: app + coaching; -1: app alone)
- **R**: response status (1: responder; 0: nonresponder)
- **nonresponse_week**: week classified as nonresponder (0: responders; 2,4,8: nonresponders)
- **Z2**: secondary intervention (1: vigorous step-up; -1: moderate step-up)
- **A**: intervention based on microrandomization (1: message sent, 0: message not sent)
- **days_nonrespond**: days since entering the MRT (vary based on nonresponse_week)
- **outcome**: binary variable representing food tracking within 12 hours of microrandomization (1: food tracking occurred; 0: food tracking did not occur)
- **prob**: probability of intervention with microrandomization (2/3: nonresponders; 0: responders)

### sim_data_for_distal_nonresponders_only.R

This file contains simulated data for the distal analysis. It includes only nonresponders. The variables included are detailed below:

- **id**: participant id
- **female**: binary sex (1: female; 0: male)
- **sex_centered**: binary sex (as defined above) mean centered
- **bmi_bl**: body mass index (bmi) at baseline in lbs/in^2
- **bmi_bl_centered**: bmi mean centered
- **meanAcentered**: average number of interventions received per participant, mean centered
- **Z1**: initial intervention (1: app + coaching; -1: app alone)
- **Z2**: secondary intervention (1: vigorous step-up; -1: moderate step-up)
- **weight_change_distal**: difference in weight from baseline to 6 months (baseline minus 6 months)

### Simulated_Data_Analysis_Final.R

This file contains code for the analysis of proximal and distal effects in a hybrid SMART-MRT design. It pulls from the file "Hanna_Estimator_Resp.R" to perform the proximal analysis and uses the two datasets "sim_data_for_proximal.R" and "sim_data_for_distal_responders_only.R".

### Estimator_For_Proximal_Analysis.R

This file contains code for conducting the proximal effect analysis with a binary proximal outcome. The code is based on Qian et al., 2021.
Qian, T., Yoo, H., Klasnja, P., Almirall, D., & Murphy, S. A. (2021). Estimating time-varying causal excursion effects in mobile health with binary outcomes. Biometrika, 108(3), 507-527.


