# Understanding wellbeing

This repository hosts the primary repository of the ffacoach, the result of the New York University MSBA Capstone project of 2020. 

### [Access](https://webapp.facharztpraxis-fuer-allgemeinmedizin.de/)

## The Project | Summary

Modern technology can substantially contribute to improving the experience for both patients and medical doctors within the health care system addressing key issues industrial nations are exposed to nowadays. The main problem we are focused on is that practitioners feel forced to treat more patients in shorter periods of time due to cost pressures; this is especially the case in countries such as Germany where a large share of the population is over 65 leading to many chronic treatments and cost increases. As a result of these pressures, general health practitioners are at times forced to adopt a practice referred to as “5-minute-medicine”, a practice that forces appointments to be shorter. If patient consultations are expedited there is increased risk of misdiagnosis, mistreatment and patient discomfort. How we aim to solve this problem is to supplement medical services through an automated Health Coach in the form of a mobile application. The primary objective of that application is to give patients the ability to quantify their wellbeing and support their self-management with actionable recommendations. Simultaneously, this app shall give doctors the option to better understand the wellbeing of their patients accessing this data when granted for more transparency and objectivity about the patients’ lifestyle. 

The development of the application is predicated on predictive, causal, and statistical analysis as well as software engineering. A representative survey of >24,000 participants and 269 variables in Germany primarily poses the data foundation for the development. It is the result of a study undertaken by the Robert-Koch-Institute (RKI) between 2014 and 2015. Upon this, amongst a diverse set of models, we conducted tests using ordinal logistic regressions and random forest to better understand drivers of wellbeing and marginal effects of simulated improvement as the foundation for recommendations. As a first important step, we could validate general patterns and drivers for human wellbeing against state of research, medical knowledge and intuition. Overall, we analyzed over 90 drivers for wellbeing across six overarching themes describing socio-demographics, health history, preventive measures taken or lifestyle attributes of participants. Analyzing and validating existing drivers on wellbeing, we conducted both propensity score matching and instrumental variable test to assess if causal relationships existed between visiting a medical practitioner and patient wellbeing. Understanding driving patterns of wellbeing, we focused on creating recommendations to leading a healthier life on a very concrete and actionable yet personalized basis. This recommender system inferred marginal effects on the expected wellbeing score, simulating improvement on actionable drivers. We concluded that from a marginal effect standpoint, as expected, a one-unit change to one of the variables (e.g., following one recommendation, ceteris paribus) will result in a small increment of improvement. Yet, changing just one variable will not result in your wellbeing score increasing materially, which is why issuance of various guiding suggestions targeted to an individual is useful. Furthermore, we attempted to transform a gain in wellbeing into expected increase in life-expectancy as a more tangible metric.

This project is supported through a partnership with a general practitioner in Germany – the Facharztpraxis fuer Allgemeinmedizin (FFA). With approx. 11,000 patients, FFA specializes in General Medicine, Sports Medicine and Emergency Care. This partnership allowed us to obtain key insights in the German health care market as well as access to medical knowledge and staff. Our partnership with FFA will continue beyond the scope of this project into the testing phase; with the guidance of FFA, their patients are invited to participate. For this and in close collaboration with FFA, we will continue to address various data engineering challenges, amongst others accounting for and complying with data privacy regulation within the EU and appropriate certification. This is crucial, specifically as it relates to both sharing approved patient information with doctors and storing such sensitive patient information securely.

In order to conclude, although additional tests are required, we believe that the effective use of the Health Coach will assist both patients and health care experts in better managing their wellbeing and their health practices respectively. Furthermore, it might even be a feasible solution to gather valuable longitudinal data for research purposes monitoring infected patients during the SARS-CoV-2 (COVID19) pandemic complementing existing solutions for infection tracing.

## The application repositories: 

* [Application:](https://github.com/justusfowl/ffa-app)
* [Backend app layer:](https://github.com/justusfowl/ffa-app-be)
* [Analytical layer (Python):](https://github.com/justusfowl/ffa-app-mlbe)

## Data repository: 

* [Google Drive](https://drive.google.com/drive/u/1/folders/1wvAqAupiXMe-EwNtcfGe2nSTMULleidR)

## Software 

The application makes use of several open source packages and is built for research, non-commercial use. 

## Contributions
* [A. Aybar](https://github.com/andresaybar)
* [J. Williams](https://github.com/williamsjerimiah)
* [J. Asmar](https://github.com/joseasmar)
* [W. Fluney](https://github.com/Fluneb)
* [U. Kaulfuss](https://github.com/justusfowl)