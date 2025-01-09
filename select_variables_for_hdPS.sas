/************************************************************************************************************************/
/* Program: /local/projects/medicare/sglt_cvs/programs/4c_hc_covariates.sas                                             */
/* Purpose:  select variables for hdPS 	                                                                                */
/* Author: Virginia Pate                                                                                                */
/************************************************************************************************************************/


options sasautos=(SASAUTOS "/local/projects/medicare/sglt_cvs/programs/macros");
%setup(full, analysis/select_variables_for_HDPS, saveLog=Y);



%LET drug1=SGLT;
%LET drug2=GLP;

/****************************************************************************************/
/* STEP 4 (hdPS). Prioritize Covariates.                                                        .    */
/****************************************************************************************/
*STEP 1: SELECT ALL POSSIBLE CANDIDATE VARIABLES BASED ON MAX NUMBER AND MINIMUM PREVALENCE. ALSO GRAB OUTCOME VARIABLES;
%macro select_vars();
   *Add code to get count of patients with condition from the prevalence;
   proc sql noprint; select count(*) into :N_cohort from out.hdCov_&drug1.v&drug2 ; quit;

	*Select top &max_num prevalent variables from each dimension;
	%macro select_vars_sub(dimension=);
		%GLOBAL list_&dimension N_&dimension;	
		proc sql noprint;
			select distinct strip(code) into :list_&dimension separated by ' '
				from out.ref_hdCov_&drug1.v&drug2 /*generated from the 4c_hd_covariates.sas*/
				where lowcase(dimension)="&dimension" and prev_order<=&max_num and prevalence>=&min_prev/100 and count>=100;
			%LET N_&dimension = &SqlObs;
		quit;

		%DO d=1 %TO &&&N_&dimension; %GLOBAL &dimension.&d; %END;
		proc sql noprint;
			select distinct strip(code) into :&dimension.1-:&dimension.&max_num
				from out.ref_hdCov_&drug1.v&drug2
				where lowcase(dimension)="&dimension" and prev_order<=&max_num and prevalence>=&min_prev/100 and count>=100;
		quit;

		%PUT &dimension: &&&N_&dimension; 
	%mend;

	%select_vars_sub(dimension=dx&dxGroup._inpt);*61;
	%select_vars_sub(dimension=dx&dxGroup._outpt);*200;
	%select_vars_sub(dimension=cpt5_inpt);*24;
	%select_vars_sub(dimension=cpt5_outpt);*200;
	%select_vars_sub(dimension=atc&atcGroup._outpt);*92;
*577 total;

	*Keep only those variables on the cohort dataset;
	data hdcov1_&drug1.v&drug2._prev&min_prev._num&max_num;
		set out.hdCov_&drug1.v&drug2; /*ALSO generated from the 4c_hd_covariates.sas*/
		keep bene_id indexdate
		  %DO d=1 %TO &&N_dx&dxGroup._inpt; 
		  		dx&dxGroup._inpt_once_&&&dx&dxGroup._inpt&d 
				dx&dxGroup._inpt_sporadic_&&&dx&dxGroup._inpt&d 
				dx&dxGroup._inpt_frequent_&&&dx&dxGroup._inpt&d %END;
		  %DO d=1 %TO &&N_dx&dxGroup._outpt; 
		  		dx&dxGroup._outpt_once_&&&dx&dxGroup._outpt&d 
				dx&dxGroup._outpt_sporadic_&&&dx&dxGroup._outpt&d 
				dx&dxGroup._outpt_frequent_&&&dx&dxGroup._outpt&d %END;
		  %DO d=1 %TO &&N_cpt5_inpt;  /*9/29/2024 Tadded one more "&"???*/
		  		cpt5_inpt_once_&&cpt5_inpt&d
				cpt5_inpt_sporadic_&&cpt5_inpt&d 
				cpt5_inpt_frequent_&&cpt5_inpt&d %END;
		  %DO d=1 %TO &&N_cpt5_outpt; /*9/29/2024 Tian added one more "&"???*/
		  		cpt5_outpt_once_&&cpt5_outpt&d
				cpt5_outpt_sporadic_&&cpt5_outpt&d 
				cpt5_outpt_frequent_&&cpt5_outpt&d %END;
		  %DO d=1 %TO &&N_atc&atcGroup._outpt; 
		  		atc&atcGroup._outpt_once_&&&atc&atcGroup._outpt&d
				atc&atcGroup._outpt_sporadic_&&&atc&atcGroup._outpt&d 
				atc&atcGroup._outpt_frequent_&&&atc&atcGroup._outpt&d%END;
			;
	run;
	/*----------------9/18/2024, according to the hdPS paper, 
A code that appeared above the 75th percentile number of times would have a "TRUE" value for all 3 recurrence variables. 
"If any of the values were equal, the variable representing the higher cutpoint was dropped." */
/*	        case when b.median = b.q3    then . else calculated cov_frequent end as final_cov_frequent,*/
/*	        case when b.median = 1       then . else calculated cov_sporadic end as final_cov_sporadic*/
	/*9/27/2024 -- moved variable drop section here -- actually drop variables rather than set them to missing*/

	%macro drop_vars(dimension=dx3_inpt);
		proc sql noprint;
			select distinct code into :sporadic_code1-:sporadic_code2000
				from out.ref_hdCov_&drug1.v&drug2 where median=1 and lowcase(dimension)="&dimension"
					and prev_order<=&max_num and prevalence>=&min_prev/100 and count>=100;
			%LET Nsporadic=&SqlObs;

			select distinct code into :frequent_code1-:frequent_code2000
				from out.ref_hdCov_&drug1.v&drug2 where median=q3 and lowcase(dimension)="&dimension"
					and prev_order<=&max_num and prevalence>=&min_prev/100 and count>=100;
			%LET Nfrequent=&SqlObs;
		quit;

		%GLOBAL drop_list_&dimension;
		   %LET drop_list_&dimension=;
		%DO i=1 %TO &Nsporadic; %LET drop_list_&dimension = &&drop_list_&dimension &dimension._sporadic_&&sporadic_code&i; %END;
		%DO i=1 %TO &Nfrequent; %LET drop_list_&dimension = &&drop_list_&dimension &dimension._frequent_&&frequent_code&i; %END;
	%mend;

	%drop_vars(dimension=dx&dxGroup._inpt);
	%drop_vars(dimension=dx&dxGroup._outpt);
	%drop_vars(dimension=cpt5_inpt);
	%drop_vars(dimension=cpt5_outpt);
	%drop_vars(dimension=atc&atcGroup._outpt);


	data hdcov2_&drug1.v&drug2._prev&min_prev._num&max_num;	
	 set hdcov1_&drug1.v&drug2._prev&min_prev._num&max_num;
		drop &&drop_list_dx&dxGroup._inpt 
			 &&drop_list_dx&dxGroup._outpt 
			 &drop_list_cpt5_inpt 
			 &drop_list_cpt5_outpt
			 &&drop_list_atc&atcGroup._outpt;
	run;

	*Add outcome variables;
	*TO DO: THIS STEP SHOULD CREATE OUR FINAL ANALYTIC DATASET - APPLY ADDITIONAL EXCLUSION CRITERIA AS NEEDED AND CREATE 2-YEAR OUTCOME VARIABLE;
 	proc sql;
	create table hdcov3_&drug1.v&drug2._p&min_prev._n&max_num._i&dxGroup._A&atcgroup. as
				select a.*, b.&outcome_date as outcome_date format=date9., c.&drug1.,
						 c.filldate2, c.endPartD, c.discontDate, c.switchAugmentDate_comparator, c.death_dt, c.abenddt, c.censorDate_itt,
						 c.excludeFlag_preFill2Initiator, c.excludeFlag_sameDayInitiator, c.excludeFlag_prevalentUser 
		    from hdcov2_&drug1.v&drug2._prev&min_prev._num&max_num as a 
				left join temp.outcomes_sgltvglp as b on a.bene_id=b.bene_id and a.indexdate=b.indexdate
				left join temp.newusers_sgltvglp as c on a.bene_id=c.bene_id and a.indexdate=c.indexdate;
    quit;

	 data    hdcov_&drug1.v&drug2._p&min_prev._n&max_num._i&dxgroup._A&atcgroup.;
	 	set hdcov3_&drug1.v&drug2._p&min_prev._n&max_num._i&dxgroup._A&atcgroup.;

		days_itt = min(outcome_date, censorDate_itt, (filldate2+2*365)) - filldate2 + 1;
		if days_itt<=0 then event_itt=.;
			else if filldate2<=outcome_date<=censorDate_itt then event_itt=1;
			else event_itt=0;
 
		censorDate_at = min(death_dt, discontDate, switchAugmentDate_comparator, (filldate2+2*365));
		days_at = min(outcome_date, censorDate_at) - filldate2 + 1;
		if days_at<=0 then event_at=.;
			else if filldate2<=outcome_date<=censorDate_at then event_at=1;
			else event_at=0;
	run;
%mend;


*STEP 2: CALCULATE BIAS FOR EACH CANDIDATE VARIABLE (i.e., Step 4 in HDPS paper);
%macro bias();	
	proc datasets lib=work nolist nodetails; delete bias; run; quit;

	proc sql noprint;
		select distinct name into :variable1-:variable5000
			from dictionary.columns where upcase(libname)='WORK' and upcase(memname)="HDCOV_&drug1.V&drug2._P&min_prev._N&max_num._I&dxgroup._A&atcgroup."
				and scan(upcase(name),1,'_') in ("DX&DXGROUP" "CPT5" "ATC&ATCGROUP");
		%LET N = &SqlObs;
	quit;
	

	%DO i=1 %TO &N;
		proc sql noprint;
			select distinct &&variable&i into :level1-:level2 from hdcov_&drug1.v&drug2._p&min_prev._n&max_num._i&dxgroup._a&atcgroup;
			%LET Nlevels=&SqlObs;
		quit;

		%IF &Nlevels=2 %THEN %DO;
			*Step 2a: Get RR of outcome by confounder (RRcd)-- 9/25/24 -- change to (1/rr), change implemented below;
			*RRCD: the relative risk for the univariate association between the binary covariate and the study outcome.;

		   	ods output RelativeRisks=_rr_(keep=statistic value);
				proc freq data=hdcov_&drug1.v&drug2._p&min_prev._n&max_num._i&dxgroup._A&atcgroup.;			  
					tables &&variable&i * event_itt / relrisk 
					cmh /*8/1/2024 Tian added this for Mantel-Haenszel procedure according to discussion w/ Til, Richie*/
					;
				run;

				*9/25/24 - use (1/rr) in place of rr;
				data _rr2_; 
					set _rr_ end=eof;
					retain rr;
					if statistic='Relative Risk (Column 2)' and value ne . then rr=(1/value);
					if eof and rr=. then rr=(1/value);
					if eof then output;
					keep rr;
				run;

			*Step 2b: Get prevalence of confounder by exposure level 
			P0(i.e.,Pc0)=the prevalence of the binary covariate within the unexposed group,
			P1(i.e.,Pc1)=the prevalence of the binary covariate within the   exposed group;
				ods output CrossTabFreqs=_prev_(where=(_type_='11' and &&variable&i=1) keep=_type_ &drug1 &&variable&i colPercent);
				proc freq data=hdcov_&drug1.v&drug2._p&min_prev._n&max_num._i&dxgroup._A&atcgroup.;
					tables &&variable&i * &drug1;
				run;
				
				data _prev2_; 
					set _prev_ end=eof; 
					retain p0 p1; 
						 if &drug1 =1 then p1=colPercent/100; 
					else if &drug1 =0 then p0=colPercent/100;
					if eof then do; 
						 if p0=. then p0=0; 
						 if p1=. then p1=0; output; 
					end; 
					keep p0 p1; 
				run;

			*Step 2c: Calculate bias as a function of RR (step 2a) and P1, P0 (step 2b);
			data _bias_;*(keep=variable log_bias); 
				length variable $50; variable="&&variable&i"; 
				merge _rr2_ _prev2_;
			/*7/30/2024 Tian edited: according to the Erratum of hdPS: https://pubmed.ncbi.nlm.nih.gov/29958191/
			No matter RRcd >=1 or <1*/			
				biasM = (p1*((1/rr)-1) + 1) / (p0*((1/rr)-1) + 1);
				log_bias = log( abs(biasM) );
			run;

			%IF %SYSFUNC(exist(bias)) %THEN %DO; proc append base=bias data=_bias_; run; %END;
			%ELSE %DO; data bias; set _bias_; run; %END;

			proc datasets lib=work nolist nodetails; delete _rr_ _rr2_ _prev_ _prev2_ _bias_; run; quit;
		%END;
	%END;
%mend;


/****************************************************************************************/
/* STEP 5 (hdPS). Select Covariates.                                                    */
/****************************************************************************************/
/*%LET dxgroup=3; %LET atcgroup=3; %LET max_num=200; %LET min_prev=1; %LET outcome_date=hfprimary_icd10dx_date;*/
%macro run_hdps(dxGroup=3, atcGroup=3, max_num=100, min_prev=1, outcome_date=HFPRIMARY_ICD10DX_date);
	*Step 1: Select candidate covariates (step 4 from HDPS paper);
	%select_vars()

	*Step 2: Calculate bias (step 4 from HDPS paper);	
	%bias()

	*Step 3: Select final variables (step 5 from HDPS paper);
	proc sort data=bias; by descending log_bias; run;

	data out.hdps_vars_&drug1.v&drug2._p&min_prev._n&max_num._i&dxgroup._A&atcgroup.;
		set bias; 
		priority_order=_N_;
	run;
	
	data top500;
    	set out.hdps_vars_&drug1.v&drug2._p&min_prev._n&max_num._i&dxgroup._A&atcgroup.;
    	if _N_ <= 500;
	run;
	
	proc sql noprint;
	    select distinct variable into :var1-:var500 from top500;
	quit;

	proc sql;
	    create table hdps0_&drug1.v&drug2._p&min_prev._n&max_num._i&dxgroup._A&atcgroup. as
	    select bene_id, indexdate, days_itt, event_itt, days_at, event_at, 
		    outcome_date, censorDate_itt, filldate2, death_dt, 
		    censorDate_at, discontDate, switchAugmentDate_comparator, 
		    excludeFlag_preFill2Initiator, 
		    excludeFlag_sameDayInitiator,
		    excludeFlag_prevalentUser, 
		    %DO i=1 %TO 500; &&var&i %IF &i<500 %THEN , ; %END;
	    from hdcov_&drug1.v&drug2._p&min_prev._n&max_num._i&dxgroup._A&atcgroup.;
	quit;

	proc sql;
	    create table ana.hdps_&drug1.v&drug2._p&min_prev._n&max_num._i&dxgroup._A&atcgroup. as
	    select a.*, b.&drug1., b.age, b.sex, b.race
	    from hdps0_&drug1.v&drug2._p&min_prev._n&max_num._i&dxgroup._A&atcgroup. 
	    as a left join temp.newusers_&drug1.v&drug2. as b 
	    on a.bene_id=b.bene_id and a.indexdate=b.indexdate;
	quit;

%mend;
			
/*primary analysis*/			
%run_hdps(dxGroup=3, atcgroup=3,  min_prev=1, max_num=200, outcome_date=HFPRIMARY_ICD10DX_date);

/*sensitivity analysis 1*/

%run_hdps(dxGroup=3, atcgroup=3,  min_prev=1, max_num=100, outcome_date=HFPRIMARY_ICD10DX_date);

/*sensitivity analysis 2*/
%run_hdps(dxGroup=3, atcgroup=3,  min_prev=2, max_num=200, outcome_date=HFPRIMARY_ICD10DX_date);

/*sensitivity analysis 3*/
%run_hdps(dxGroup=3, atcgroup=3,  min_prev=5, max_num=200, outcome_date=HFPRIMARY_ICD10DX_date);

/*sensitivity analysis 4*/
%run_hdps(dxGroup=4, atcgroup=3,  min_prev=1, max_num=200, outcome_date=HFPRIMARY_ICD10DX_date);

/*sensitivity analysis 5*/
%run_hdps(dxGroup=3, atcgroup=4,  min_prev=1, max_num=200, outcome_date=HFPRIMARY_ICD10DX_date);



/*skip the %bias macro, i.e. no variable selection for top 500, generate a dataset*/

/*%run_hdps(dxGroup=3, atcgroup=4,  min_prev=1, max_num=200, outcome_date=HFPRIMARY_ICD10DX_date);*/

/*7/31/2024: step 5 from HDPS continue, select top 500*
data top500;
    set out.hdps_vars_sgltvglp_p1_n200_i3_a3;
    if _N_ <= 500;
run;*/

/*BE CAREFUL!!!: Need to clean work libarary before running %run_hdps due to this step:
data bias; set bias:; run; which keep stacking datasets leading to duplicates!!!
proc sql; select count(distinct variable) as distinct_var from top500; quit;
proc sql; select count(         variable) as          var from top500; quit;
proc sql; select count(distinct variable) as distinct_var from ana.hdps_vars_sgltvglp_p1_n200_i3_a3;run;
proc sql; select count(		    variable) as 		  var from ana.hdps_vars_sgltvglp_p1_n200_i3_a3;run;
proc contents data=out.hdCov_&drug1.v&drug2;run;*/
/*primary analysis & sensitivity anlaysis 1 & 2*/
/*%run_hdps(dxGroup=3, atcgroup=3,  min_prev=1, max_num=200, outcome_date=HFPRIMARY_ICD10DX_date);*/
/*/*Sensitivity analysis 3: First, we utilized the 100 most prevalent codes (n=100) in each dimension.*/*/
/*%run_hdps(dxGroup=3, atcgroup=3,  min_prev=1, max_num=100, outcome_date=HFPRIMARY_ICD10DX_date) */
/*/*Sensitivity analysis 4: Second, we applied a cutoff c=0.02 for these codes.*/ */
/*%run_hdps(dxGroup=3, atcgroup=3,  min_prev=2, max_num=200, outcome_date=HFPRIMARY_ICD10DX_date);*/
/*/*Sensitivity analysis 5: Second, we applied a cutoff c=0.05 for these codes.*/ */
/*%run_hdps(dxGroup=3, atcgroup=3,  min_prev=5, max_num=200, outcome_date=HFPRIMARY_ICD10DX_date);*/
/*/*Sensitivity analysis 6: Third, we employed a 4-digit convention for the ICD-10 codes.*/*/
/*%run_hdps(dxGroup=4, atcgroup=3,  min_prev=1, max_num=200, outcome_date=HFPRIMARY_ICD10DX_date); */
/*/*Sensitivity analysis 7: Fourth, we used the 4th level for the ATC codes. */*/
/*%run_hdps(dxGroup=3, atcgroup=4, min_prev=1, max_num=200, outcome_date=HFPRIMARY_ICD10DX_date);*/
/*/*Sensitivity analysis 8: Fourth, we used the 4th level for the ATC codes. */*/
/*%run_hdps(dxGroup=3, atcgroup=4, min_prev=2, max_num=200, outcome_date=HFPRIMARY_ICD10DX_date);*/
/*/*Sensitivity analysis 9: Fourth, we used the 4th level for the ATC codes. */*/
/*%run_hdps(dxGroup=3, atcgroup=4, min_prev=5, max_num=200, outcome_date=HFPRIMARY_ICD10DX_date);*/
/*/*Sensitivity analysis 10: Fourth, we used the 4th level for the ATC codes. */*/
/*%run_hdps(dxGroup=3, atcgroup=4, min_prev=5, max_num=100, outcome_date=HFPRIMARY_ICD10DX_date);*/
/*/*Sensitivity analysis 11: Fourth, we used the 4th level for the ATC codes. */*/
/*%run_hdps(dxGroup=3, atcgroup=4, min_prev=5, max_num=50, outcome_date=HFPRIMARY_ICD10DX_date);*/
