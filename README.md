# webcovid19
Suplementar material for manuscript entitled *Estimating the Brazilian states’ demand for intensive care unit and clinical hospital beds in COVID-19 pandemic: prediction model development*, São Paulo Medical Journal. Reference **SPMJ-2020-0517**.


Consider a compartmental model21,22 in which the population is divided into susceptible (S), exposed (E), infectious (I), and recovered (R) individuals. model dynamics is illustrated in Figure 1.

![Figure 1](https://github.com/joaoflavioufmg/webcovid19/blob/main/fig1.png)

A rate of β of S individuals in contact with I (S-I contact) becomes E and progress in the course of the incubation period in a rate σ to state I. While a rate γ of I recover from the disease, a rate of I evolves to death. A fraction ξ of R recovered individuals may become (or not, if ξ = 0) re-susceptible (S), therefore, an SEIRS model. The effect of testing a segment of the population is modeled by introducing the rate of transmission for individuals with detected infections (β<sub>D</sub>), detected exposed (D<sub>E</sub>), and detected infectious (D<sub>I</sub>) states. Let ψ<sub>E</sub> and ψ<sub>I</sub> be the probability of positive tests for exposed and infected individuals, and Q, the rate of individuals with detected infection interacting with the population. Then D<sub>E</sub> and D<sub>I</sub> results from θ<sub>E</sub>ψ<sub>E</sub> and θ<sub>I</sub>ψ<sub>I</sub> rates of testing exposed E and infected I individuals, respectively. The model describes the full spectrum of disease. Let N be the estimate of an affected population, thus, N=S+E+I+D_E+D_I+R. Let F be the number of infected-related fatalities. The model states are obtained by equations (1)-(7).

![Model's equations (1)-(7)](https://github.com/joaoflavioufmg/webcovid19/blob/main/eq-seirs.png)

To represent the ICU and clinical beds’ dynamics, consider that a fraction of infected individuals is asymptomatic (I<sub>A</sub>). The rate α for symptomatic cases (I<sub>S</sub>=I-I<sub>A</sub>) requiring hospitalization (H) is \alpha I<sub>S</sub>. Let T be the number of planning days during the pandemic with t corresponding to each admission day on the hospital. For clinical bed dynamics, consider L<sub>s</sub> and L<sub>f</sub> as the average length of stay (LoS) of surviving patients,  and the patients that passed away, respectively. For surviving patients on ICU beds, the average length of stay is L<sub>b</sub>+L<sub>d</sub>+L<sub>a</sub>, where L<sub>b</sub> represents the surviving patients of clinical beds re-directed to ICU beds, L<sub>d</sub> is the average LoS of surviving patients in ICU beds, and L<sub>a</sub>, as the average LoS in clinical beds of surviving patients after being re-directed from ICU beds. For deceased victims on ICU bed dynamics, the average LoS is L<sub>c</sub>+L<sub>i</sub>, representing the average period a deceased patient will stay in clinical and ICU beds, respectively. Given the average length of stay metrics for ICU and clinical dynamics, we calculate admission and leave days for each patient profile. For a surviving clinical patient, whose admission day is t, the expected day that the patient leaves the clinical bed is T<sub>o</sub>=t+L<sub>s</sub>–1. A clinical patient that progressed to death is expected to be removed from the hospital in T<sub>f</sub>=t+L<sub>f</sub>-1. A surviving ICU patient is admitted in T<sub>i</sub>=t+L<sub>b</sub>. In T<sub>r</sub>=T<sub>i</sub>+L<sub>d</sub>, the patient return to the clinical bed, and in T<sub>c</sub>=T<sub>r</sub>+L<sub>a</sub>-1, the patient leaves the clinical bed. For a deceased ICU patient, the admission day is T<sub>d</sub>=t+L<sub>c</sub> , and T<sub>u</sub>=T<sub>d</sub>+L<sub>i</sub>-1 is the expected day that the patient is removed from the hospital. Let H<sub>t</sub> be the daily admission of patients in hospitals, and τ the fraction of hospitalized cases that require critical care in ICU. Also, consider ζ and η the rates of critical patients that evolve to death on ICU beds and clinical beds, respectively. Therefore, the daily use of ICU beds and clinical beds, represented by U<sub>tICU</sub> and U<sub>tC</sub>, is obtained by the equations (8).

![Beds dynamics equations](https://github.com/joaoflavioufmg/webcovid19/blob/main/eq-beds.png)

According to equation (8) the ICU bed patients follow two paths, survival and decease. Surviving patients rest in ICU beds from the day they enter the hospital until the day they return to clinical beds. Deceased patients rest in ICU beds from the day they are transferred from the clinic (clinical-ICU) until the day the deceased is removed from hospital. Equation (9) corresponds to the use of clinical beds from patients following survival and decease paths. For the dynamics of the clinical path, surviving patients use clinical beds from the day they enter to hospital until the expected day the patient leaves the clinical bed. Deceased patients also use clinical beds from the day they enter to hospital until the expected day the deceased is removed from the hospital. For the ICU paths dynamics, surviving ICU patients use clinical beds from the day they enter -to hospital until the admission day in ICU, and after ICU stay. Surviving ICU patients use clinical beds from the day the ICU patient returns to clinical bed until the expected day the patient leaves the hospital. Finally, deceased ICU patients use clinical beds from the day they enter to hospital until the admission day in ICU. For a given Brazilian state, the available ICU beds (A<sub>ICU</sub>) and available clinical beds (A<sub>C</sub>) for COVID-19’s infected patients correspond to a fraction of a give states’ clinical and ICU beds. The daily availability of ICU and clinical beds for a given time t are described in equations (10) and (11). 

Figures 2 and 3 show that the model is periodically adjusted to the real number of cases providing a fair estimate of future infections. Such information feeds the dynamic beds equations that forecast the reduction of ICU and ordinary beds capacity of the state health system. In May 2020, the ICU beds’ capacity for this state is different from August 2020. Thus, we adopted a simplified assumption considering a single capacity increase on July/2020.

![Figure 2](https://github.com/joaoflavioufmg/webcovid19/blob/main/fig3.png)

![Figure 3](https://github.com/joaoflavioufmg/webcovid19/blob/main/fig2.png)
