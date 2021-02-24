{smcl}
{* 17feb2021} {...}

{hline}
help for {hi:sumindex}
{hline}

{title: Summary Index}

{p 4 8 2}{cmd:sumindex} [{it:varlist}] {cmd:,}
	{cmdab:gen:erate(}{it:varname}{cmd:)} [{cmdab:b:ase(}{it:condition}{cmd:)}
	{cmdab:r:eplace} {cmdab:norm:alise} {cmdab:nop:airwise}] {p_end}

{p 4 4 2}
{cmd:by} {it:...} {cmd::} may be used with {cmd:sumindex}; see {help by}.

{title:Description}

{p 4 4 2}
{cmd:sumindex} generates a summary index variable called {it:varname} using {it:varlist} as the index components. The procedure implements the summary index described by Anderson (2008) in Section 3.2.1 and more formally in Appendix A, with several additional modifications. This summary index is a weighted mean of several standardised outcomes where the weights come from the inverse of the covariance matrix. The standardisation uses a base (or reference) group which is supplied in the {cmdab:b:ase(}{it:condition}{cmd:)} option.

{p 4 4 2}
The two required conditions for the variables in {it:varlist} are that they are numerical and that the values are non-missing for at least two observations (which is required to calculate the standard deviation). If the variables supplied in {it:varlist} fail any of these conditions, a missing value for the index will be returned. In case only one index component is supplied, this code will simply create a standardised variable of the singular index component. Note that the variables should be created such that higher values indicate a "better" outcome.

{p 4 4 2}
The base group identifies the group of observations that act as the reference group for normalisation. Anderson (2008) normalises the index component variables by demeaning the variables and then dividing by the standard deviation of the base group. Furthermore, the covariance matrix is calculated using the entire sample. If the option {cmdab:norma:lise} is selected, this will modify the index in three ways. Firstly, instead of demeaning each index component variable, the {cmdab:norma:lise} option will subtract the base group variable mean from each index component variable. Secondly, the covariance matrix is determined using only observations from the base group. Thirdly, after demeaning using the base group mean, if there are missing observations the final index will not have a mean of zero and a standard deviation of one for the base group. In order to address this, the {cmdab:norma:lise} option will normalise the index by subtracting the base group index mean and dividing by the base group standard deviation.

{p 4 4 2}
The {cmdab:nop:airwise} option allows for the calculation of the covariance matrix using only the set of observations for which there are no missing values for all index component variables. Without the {cmdab:nop:airwise} option specified, each pairwise covariance in the covariance matrix is calculated using the sample of non-missing observations for the pair of variables.

{title:Options}

{p 4 4 2}
{cmdab:gen:erate(}{it:varname}{cmd:)} specifies the name of the new index variable that will be created.

{p 4 4 2}
{cmdab:b:ase(}{it:condition}{cmd:) specifies the condition to identify the base group. If this option is not supplied, the entire sample is considered to be the base group.

{p 4 4 2}
{cmdab:r:eplace} specifies that if the variable {it:varname} already exists, then {cmd:sumindex} will replace the values in the variable.

{p 4 4 2}
{cmdab:norm:alise} specifies that the index is calculated i) using the base group mean to normalise the index component variables, ii) using a covariance matrix that is calculated using observations from the base group and iii) is normalised such that the index will have a mean of zero and a standard deviation of one for the base group.

{p 4 4 2}
{cmdab:nop:airwise} specifies that the covariance matrix is created using only the set of observations for which there are no missing values for all index component variables.

{title:Examples}

{p 4 4 2}
Suppose you want to make an economic index ({it:econ_index}) with the following index components: {it:income}, {it:assets}, {it:consumption}. The base group is defined as the control group identified with a dummy variable called {it:control}. Additionally, the data has multiple rounds identified with the variable {it:round}. Therefore, the following code would be run:

{cmd:. bysort round: sumindex income assets consumption, gen(econ_index) base(control == 1)}

{title:Author}

    {bf:Paul Brimble}
    Blavatnik School of Government
    University of Oxford
    paul.brimble@bsg.ox.ac.uk
    http://pbrimble.github.io

{title:References}

{p 4 4 2}
Michael L. Anderson (2008) Multiple Inference and Gender Differences in the Effects of Early Intervention: A Reevaluation of the Abecedarian, Perry Preschool, and Early Training Projects, Journal of the American Statistical Association, 103:484, 1481-1495, DOI: 10.1198/016214508000000841
