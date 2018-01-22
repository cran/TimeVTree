\name{alcohol}
\alias{alcohol}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Alcohol Consumption Data
}
\description{
This data set contains subjects' age at the time of death, and alcohol drinking habits.
The data set includes 7990 subjects and 7610 events.
}
\usage{
data('alcohol')
}
\format{
\itemize{
  \item{time}{:
Subject's age at death (possibly right censored)
}
  \item{event}{:
Outcome indicator.\cr 1 = death\cr 0 = censored\cr
}
  \item{alc}{\cr
0 =  no alcohol consumption\cr
1 = moderate alcohol consumption\cr
4 = excessive alcohol consumption
}
}
}
\source{
Data is from the Honolulu Heart Program/Honolulu Asia Aging Study (HHP/HAAS).
The HHP/HAAS was reviewed and approved by the Kuakini Hospital IRB, Kuakini
    Hospital, Honolulu, HI.
}
