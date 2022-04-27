/* SAS macros for performing 0-1 augmented, 0-augmented, 1-augmented, and non-augmented lemon-squeezer beta regression models. */

/* The *Beta_Regression* macro is the highest level macro and is where the user could define the model desired, initial starting values, and the optimization technique. */ 
/* This is the only macro that required input by the user. */


*Call the highest level macro;
*Beta_Regression(Dataset,tech,details,mu_vars,phi_vars,zero_vars,one_vars,depvar,params,bounds);
*;
*dataset=the data being used;
*tech=NLMIXED optimization technique;
*details=any other options wanted on the Proc NLMIXED line;
*mu_vars=the variables we wish to use to model the change in the mean;
*phi_vars=the variables we wish to use to model the changes in phi (use 1 if not using phi);
*pizero_vars=the variables we wish to use to predict response of zero;
*pione_vars=the variables we wish to use to predict response of one;
*depvar=the dependent variable;
*params=initial starting values (WARNING:SAS default is to set initial values to 1. Choose carefully so as to not exceed the defined or natural bounds of the problem)
*Bounds=add bounds inherent in the model such as variance>0;

%Beta_Regression(Work.Propdata, trureg,
,
gender age hba1c_d smoker Tooth2 Tooth3 Tooth4,
1,
gender age hba1c_d smoker Tooth2 Tooth3 Tooth4,
gender age hba1c_d smoker Tooth2 Tooth3 Tooth4,
Prop, 
b0 0 b1 0 b2 0 b3 0 b4 0 b5 0 b6 0 b7 0 
/*d0 0 d1 0 d2 0 d3 0 d4 0 d5 0  d6 0 d7 0 */  /*if using phi, uncomment d0-d7*/
zero0 -.1 zero1 0 zero2 0 zero3 0 zero4 0 zero5 0 zero6 0 zero7 0 
one0 -.1 one1 0 one2 0 one3 0 one4 0 one5 0 one6 0 one7 0  
var phi .1,
var phi>0);

/* The SAS macro below determines which of the four models to use based on the number of parameters defined. This macro is a second level macro, 
/* and passes the initially defined parameters into the desired model. 

/* Macros to assign whether to use 0-1 augmented, 0 augmented, 1 augmented, or non-augmented model */

%Macro Beta_Regression(Dataset,tech,details,mu_vars,phi_vars,zero_vars,one_vars,depvar,params,bounds);
%if &zero_vars ne and &one_vars ne %then
%Beta_Regression_Zero_One(&Dataset,&tech,&details,&mu_vars,&phi_vars,&zero_vars,&one_vars,&depvar,&bounds,&params); *Call the ZOAB if all parameters are defined in the Beta_Regression macro;
%else %if &zero_vars ne %then
%Beta_Regression_Zero(&Dataset,&tech,&details,&mu_vars,&phi_vars,&zero_vars,&depvar,&bounds,&params); *Call the 0 augmented beta regression if pi1 is not defined in the Beta_Regression macro;
%else %if &one_vars ne %then
%Beta_Regression_One(&Dataset,&tech,&details,&mu_vars,&phi_vars,&one_vars,&depvar,&bounds,&params); *Call the 1 augmented beta regression if pi0 is not defined in the Beta_Regression macro;
%else
%Beta_Regression_Only(&Dataset,&tech,&details,&mu_vars,&phi_vars,&depvar,&bounds,&params); *Call the non-augmented model if pi0 and pi1 are not defined in the Beta_Regression macro;
%mend;

/* The following four macros performed the four regressions (0-1 augmented, 0 augmented, 1 augmented or non-augmented). */

/* Macro for the non-augmented beta regression */


%Macro Beta_Regression_Only(Dataset,tech,details,mu_vars,phi_vars,depvar,bounds, params);
%Preprocessing(&mu_vars,'b0',xb,b);*create labels for the parameters;
*%Preprocessing(&phi_vars,'d0',wd,d);
proc nlmixed data = &Dataset tech = &tech &details gconv=0;
bounds &bounds;	*specify the bounds;
parms &params; *specify the inital starting values;
mu = exp(&xb+u)/(1 + exp(&xb+u)); *calculate mu using the inverse logit;
*phi = exp(&wd);	*calculate phi using the inverse log (exp());
w = mu*phi;	*w and t are parameters used to simplify the log likelihood;
t = phi - mu*phi;
ll = lgamma(w+t)-lgamma(w)-lgamma(t)+((w-1)*log(&depvar))+((t-1)*log(1-&depvar)); *calculate the log likelihood of the beta distribution;
model &depvar ~ general(ll); *model the log likelihood;
random u ~ normal(0,var) subject=ID; *the random effect due to repeated measurements;
predict mu out=mu_results (keep=&depvar pred); *create prediction datasets for the parameters of interest;
predict phi out=phi_results (keep=&depvar pred);
run; *end of proc nlmixed;
%Postprocessing(mu_hat,mu_results,&depvar);*add labels to mu variables;
%Postprocessing(phi_hat,phi_results,&depvar); *add labels to phi variables;
data prediction2; *combine predictions into one dataset;
merge mu_results phi_results;
by record;
run;
%mend;


/* Macro for the 0 augmented beta regression */


%Macro Beta_Regression_Zero(Dataset,tech,details,mu_vars,phi_vars,zero_vars,depvar,bounds,params);
%Preprocessing(&mu_vars,'b0',xb,b);	*create labels for the parameters;
%Preprocessing(&phi_vars,'d0',wd,d);
%Preprocessing(&zero_vars,'zero0',zeroxb,zero);
proc nlmixed data = &Dataset tech = &tech &details;
bounds &bounds;	*specify the bounds;
parms &params; *specify the inital starting values;
pizero = exp(&zeroxb)/(1 + exp(&zeroxb)); *calculate pizero using the inverse logit;
mu = exp(&xb+u)/(1 + exp(&xb+u)); *calculate mu using the inverse logit;
phi = exp(&wd);	*calculate phi using the inverse log (exp());
w = mu*phi;	*w and t are parameters used to simplify the log likelihood;
t = phi - mu*phi;
if (&depvar = 0) then
ll = log(pizero);
else ll = lgamma(w+t) - lgamma(w) - lgamma(t) + ((w-1)*log(&depvar)) +
((t-1)*log(1 - &depvar)) + log(1-pizero); *calculate the log likelihood of the 0 inflated beta distribution;
model &depvar ~ general(ll); *model the log likelihood;
random u ~ normal(0,var) subject=ID; *the random effect due to repeated measurements;
predict mu out=mu_results (keep=&depvar pred); *create prediction datasets for the parameters of interest;
predict phi out=phi_results (keep=&depvar pred);
predict pizero out=pizero_results (keep=&depvar pred);
run;
%Postprocessing(mu_hat,mu_results,&depvar); *add labels to mu variables;
%Postprocessing(phi_hat,phi_results,&depvar); *add labels to phi variables;
%Postprocessing(pizero_hat,pizero_results,&depvar); *add labels to pi0 variables;
data prediction; *combine predictions into one dataset;
merge mu_results phi_results pizero_results;
by record;
run;
%mend;



/* Macro for the 1 augmented beta regression */

%Macro Beta_Regression_One(Dataset,tech,details,mu_vars,phi_vars,one_vars,depvar,bounds,params);
%Preprocessing(&mu_vars,'b0',xb,b);	*create labels for the parameters;
%Preprocessing(&phi_vars,'d0',wd,d);
%Preprocessing(&one_vars,'one0',onexb,one);
proc nlmixed data = &Dataset tech = &tech &details;
bounds &bounds;	*specify the bounds;
parms &params; *specify the inital starting values;
pione = exp(&onexb)/(1 + exp(&onexb)); *calculate pione using the inverse logit;
mu = exp(&xb+u)/(1 + exp(&xb+u)); *calculate mu using the inverse logit;
phi = exp(&wd);	*calculate phi using the inverse log (exp());
w = mu*phi;	*w and t are parameters used to simplify the log likelihood;
t = phi - mu*phi;
if (&depvar = 1) then
ll = log(pione);
else ll = lgamma(w+t) - lgamma(w) - lgamma(t) + ((w-1)*log(&depvar)) +
((t-1)*log(1 - &depvar)) + log(1-pione); *calculate the log likelihood of the 1 inflated beta distribution;
model &depvar ~ general(ll); *model the log likelihood;
random u ~ normal(0,var) subject=ID; *the random effect due to repeated measurements;
predict mu out=mu_results (keep=&depvar pred); *create prediction datasets for the parameters of interest;
predict phi out=phi_results (keep=&depvar pred);
predict pione out=pione_results (keep=&depvar pred);
run;
%Postprocessing(mu_hat,mu_results,&depvar); *add labels to mu variables;
%Postprocessing(phi_hat,phi_results,&depvar); *add labels to phi variables;
%Postprocessing(pione_hat,pione_results,&depvar); *add labels to pi1 variables;
data prediction; *combine predictions into one dataset;
merge mu_results phi_results pione_results;
by record;
run;
%mend;



/* Macro for the 0-1 augmented beta regression */

%Macro Beta_Regression_Zero_One(Dataset,tech,details,mu_vars,phi_vars,zero_vars,one_vars,depvar,bounds,params);
%Preprocessing(&mu_vars,'b0',xb,b);	*create labels for the parameters;
%Preprocessing(&phi_vars,'d0',wd,d);
%Preprocessing(&zero_vars,'zero0',zeroxb,zero);
%Preprocessing(&one_vars,'one0',onexb,one);
proc nlmixed data = &Dataset tech = &tech &details gconv=0 maxiter=500;
bounds &bounds;	*specify the bounds;
parms &params; *specify the inital starting values;
pizero = exp(&zeroxb)/(1 + exp(&zeroxb)); *calculate pizero using the inverse logit;
pione = exp(&onexb)/(1 + exp(&onexb)); *calculate pione using the inverse logit;
mu = exp(&xb + u)/(1 + exp(&xb+u)); *calculate mu using the inverse logit;
phi = exp(&wd);	*calculate phi using the inverse log (exp());
w = mu*phi;	*w and t are parameters used to simplify the log likelihood;
t = phi - mu*phi;
if &depvar = 0 then
ll = log(pizero);
else if &depvar = 1 then
ll = log(pione);
else ll = lgamma(w+t) - lgamma(w) - lgamma(t) + ((w-1)*log(&depvar)) +
((t-1)*log(1 - &depvar)) + log(1-pizero-pione);	*calculate the log likelihood of the 0-1 inflated beta distribution;
model &depvar ~ general(ll); *model the log likelihood;
random u ~ normal(0,var) subject=ID; *the random effect due to repeated measurements;
predict mu out=mu_results (keep=&depvar pred); *create prediction datasets for the parameters of interest;
predict phi out=phi_results (keep=&depvar pred);
predict pizero out=pizero_results (keep=&depvar pred);
predict pione out=pione_results (keep=&depvar pred);
run;
%Postprocessing(mu_hat,mu_results,&depvar); *add labels to mu variables;
%Postprocessing(phi_hat,phi_results,&depvar); *add labels to phi variables;
%Postprocessing(pizero_hat,pizero_results,&depvar); *add labels to pi0 variables;
%Postprocessing(pione_hat,pione_results,&depvar); *add labels to pi1 variables;
data prediction; *combine predictions into one dataset;
merge mu_results phi_results pizero_results pione_results;
by record;
run;
%mend;


/* The Preprocessing macro added labels to the independent variables in the form b0,b1,etc. */

%Macro Preprocessing(Vars,b0,xb2,b);
data HPG;
%global &xb2;
length &xb2 $200.;
&xb2=&b0;
%if &Vars ne '' %then %do;
%let n=1;
var&n="%scan(&Vars,&n,' ')";
%do %while( %scan(&Vars,&n,' ') ne );
%let n=%eval(&n+1);
var&n="%scan(&Vars,&n,' ')";
%end;
%let n_1=%eval(&n-1);
array xbv {*} $ var1--var&n_1;
%do j=1 %to &n_1;
&b&j= "&b&j";
%end;
%let one=1;
array &b{*} $8 &b&one--&b&n_1 ;
array p{1} $ 8 ('+');
array m{1} $ 8 ('*');
do i=1 to dim(xbv) while (xbv{i} ne '');
&xb2= cats(of &xb2 p{1} &b{i} m{1} xbv{i});
end;
%end;
call symput("&xb2",&xb2);
run;
%mend;


/* The Postprocessing macro inputs in the predictions, and then output the predictions with new variable names. This was used to combine the datasets within the model macros. */

%Macro Postprocessing(hat, predict, depvar);
data &predict;
retain record &depvar &hat;
set &predict;
if &depvar = . then &hat = .;
else &hat = pred;
record=_n_;
keep record &depvar &hat;
run;
%mend;



