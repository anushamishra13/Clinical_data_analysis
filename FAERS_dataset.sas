/* Step 1: Import FAERS Dataset */

PROC IMPORT DATAFILE="/home/u64494097/sasuser.v94/fda_adverse_events_2015_2026_CLEAN.csv"
    OUT=work.faers_raw
    DBMS=CSV
    REPLACE;
    GETNAMES=YES;
RUN;


/* Step 2: Clean Data */

DATA work.faers_clean;
    SET work.faers_raw;

    IF suspect_drug = "" THEN DELETE;
    IF primary_reaction = "" THEN DELETE;

    suspect_drug = UPCASE(suspect_drug);
    primary_reaction = UPCASE(primary_reaction);
RUN;


/* Step 3: Create AE Dataset */

DATA work.ae_dataset;
    SET work.faers_clean;

    KEEP report_id suspect_drug primary_reaction serious 
         is_fatal is_hospitalized is_life_threat 
         patient_age_years patient_sex country;
RUN;


/* Turn off large output */
ODS RESULTS OFF;


/* Step 4: Reaction Frequency (FIXED SYNTAX) */

PROC FREQ DATA=work.ae_dataset NOPRINT;
    TABLES primary_reaction / OUT=reaction_freq;
RUN;


/* Step 5: Serious vs Non-serious */

PROC FREQ DATA=work.ae_dataset NOPRINT;
    TABLES serious*primary_reaction / OUT=serious_freq;
RUN;


/* Step 6: Fatal cases */

PROC FREQ DATA=work.ae_dataset NOPRINT;
    WHERE is_fatal = 1;
    TABLES suspect_drug / OUT=fatal_drug_freq;
RUN;


/* Step 7: Drug-event counts */

PROC SQL;
    CREATE TABLE drug_event_counts AS
    SELECT 
        suspect_drug,
        primary_reaction,
        COUNT(*) AS count
    FROM work.ae_dataset
    GROUP BY suspect_drug, primary_reaction
    HAVING count > 50;
QUIT;


/* Turn output back ON */
ODS RESULTS ON;


/* Step 8: Display top results */

PROC SORT DATA=drug_event_counts;
    BY DESCENDING count;
RUN;

PROC PRINT DATA=drug_event_counts (OBS=20);
RUN;


/* ========================= */
/* PRR CALCULATION (FIXED)   */
/* ========================= */


/* Total reports */

PROC SQL;
    CREATE TABLE total_counts AS
    SELECT COUNT(*) AS total 
    FROM work.ae_dataset;
QUIT;


/* a = drug-event count */

PROC SQL;
    CREATE TABLE a_counts AS
    SELECT 
        suspect_drug,
        primary_reaction,
        COUNT(*) AS a
    FROM work.ae_dataset
    GROUP BY suspect_drug, primary_reaction
    HAVING a > 100;
QUIT;


/* total per drug */

PROC SQL;
    CREATE TABLE drug_totals AS
    SELECT 
        suspect_drug,
        COUNT(*) AS total_drug
    FROM work.ae_dataset
    GROUP BY suspect_drug;
QUIT;


/* total per reaction */

PROC SQL;
    CREATE TABLE reaction_totals AS
    SELECT 
        primary_reaction,
        COUNT(*) AS total_reaction
    FROM work.ae_dataset
    GROUP BY primary_reaction;
QUIT;


/* Final PRR calculation */

PROC SQL;
    CREATE TABLE prr_calc AS
    SELECT 
        a.suspect_drug,
        a.primary_reaction,
        a.a,
        d.total_drug,
        r.total_reaction,
        t.total,

        (a.a / d.total_drug) / 
        ((r.total_reaction - a.a) / (t.total - d.total_drug)) AS PRR

    FROM a_counts a
    JOIN drug_totals d 
        ON a.suspect_drug = d.suspect_drug
    JOIN reaction_totals r 
        ON a.primary_reaction = r.primary_reaction
    CROSS JOIN total_counts t;
QUIT;
PROC SQL;
SELECT COUNT(*) AS total_records FROM work.ae_dataset;
QUIT;


/* Show top PRR signals */

PROC SORT DATA=prr_calc;
    BY DESCENDING PRR;
RUN;

PROC PRINT DATA=prr_calc (OBS=20);
RUN;

/* this is for downloading the the dataset in csv */
PROC EXPORT DATA=drug_event_counts
    OUTFILE="/home/u64494097/sasuser.v94/drug_event_counts.csv"
    DBMS=CSV
    REPLACE;
RUN;
