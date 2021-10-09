{smcl}
{* 17feb2021} {...}

{hline}
help for {hi:pairmatch}
{hline}

{title: Pairwise Nearest Matching }


{p 4 8 2}{cmd:pairmatch} [{it:varlist}] {cmd:using} {it:filename} {cmd:,}
	{cmdab:var_near(}{it:varname}{cmd:)} [ {cmdab:lim:it(}{it:real}{cmd:)}
	{cmdab:sign:ed(}{it:string}{cmd:)} {cmdab:order:ing(}{it:string}{cmd:)}
	{cmdab:g:enmatch(}{it:newvarname}{cmd:)} {it:mergeoptions}] {p_end}


{title:Description}

{p}{cmd:pairmatch} compares all pairwise combinations of observations in the master and using datasets and identifies nearest pairs using the numeric variable {it:var_near}.  {cmd:pairmatch} was designed as a way to merge
tables that have rounded or approximated values on the merging variable of interest and requires a unique pairwise match. This means that an observation in the master or using dataset cannot be paired with more than one other observation. {p_end}

{p}The observations in the master dataset are matched with observations in the using dataset with the smallest absolute difference between each {it:var_near} value.{{p_end}}

{p} The {it:var_near} must be a numeric variable and variables may be specified in an optional {it:varlist} that are treated as standard merge variables in order to have pairwise matching within subsets defined by the varlist.

{title:Options}

{p 0 4}{cmd:var_near()} is required and specifies the variable in the master and using datasets that is to be
matched as closely as possible. {it:varlist} {cmd:var_near()} must be unique in both the master and using dataset.{p_end}

{p 0 4}{cmd:limit()} is optional and specifies a limit on the maximum distance between any pair of observations in units of var_near. If the ordering option is not chosen, the default ordering is weak. {p_end}

{p 0 4}{cmd:signed()} is optional and specifies that the maximum distance between any pair of observations in units of var_near be signed. This option can take values of either {it:positive} or {it:negative}.
Positive implies that differences between the main and using data be greater than the limit while negative implies that the differences be less than the limit. If no limit is supplied and this option is chosen,
the default limit is 0. {p_end}

{p 0 4}{cmd:ordering()} is optional and specifies whether any differences with the limit be {it:weak} or {it:positive}. {p_end}

{p 0 4}{cmd:genmatch()} is optional and specifies that a new variable should be created in the master datset that identifies the
specific value of {cmd:var_near} in the using dataset that was matched.  {p_end}

{p 0 4}{cmd:mergeoptions} allows the user to specify any of the standard Stata {help merge} options (such as
{cmd:update} or {cmd: replace}).  See {help merge} for more on these options.{p_end}

{title:Algorithm}
{p 0 4}
This code achieves the following algorithm. For each group (identified with
varlist), execute the following:

  1) Calculate the differences between var_near for all pairwise
     combinations of observations from the main and using data and denote
	 these with {Ax, By} where Ax denotes observation x from the main
	 data and By denotes observation y from the using data.

  2) Identify the minimum number of distinct observations in the main and
	 using data and denote this with N. If there are equal observations
	 in both the main and using data, there will be no unmatched data.

  3) For each integer n in {1,...,N}":

	 3.1) Identify the pairwise combination with the n^th smallest
	      difference and label the elements of this combination {An, Bn}.
	 3.2) Omit all other pairwise combinations that contain either An or
	      Bn. Specifically omit all pairwise combinations {An, 'Bn} and
	 	 {'An, Bn} where ' denotes "not".

  4) For each group, you are left with a set of unique nearest pairwise
	 combinations {A1,B1}, ... ,{An,Bn}. {p_end}

{p 0 4} A note on tiebreaks: in cases where there is a tie for minimum differences, this
is only problematic if there is overlap with at least one of the elements. For
example, suppose {Ax,Bx} and {Ay,By} have the same difference, then the order
of the algorithm does not make a difference. However, if we have the same
difference for {Ax,Bx} and {Ax,By}, then the arbitrary ordering from the
'stable' option will make a difference. This is because if {Ax, Bx} is chosen,
the eventual pairing of By will be affected. The intuition is further
illustrated if we have N = 1 with one main observation and two using
observations. If the main observation is equally close to the using observation,
then an arbitrary decision is made. {p_end}

{title:Example}


{title:Authors}
{p 2 6 4}Current version of {bf:pairmatch} is written and maintained by:{p_end}

	{bf:Paul Brimble}
	Blavatnik School of Government
	University of Oxford
	paul.brimble@bsg.ox.ac.uk
	http://pbrimble.github.io


{title:Also See}

{p 0 19}On-line:  help for {help merge}{p_end}
