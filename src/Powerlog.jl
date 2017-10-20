module Powerlog

using Distributions, Plots

export powerlog

# !* version 1.0 -- 12/1/06 -- pbe
# !* version .9 -- 9/26/06 -- pbe
# program define powerlog, rclass
# version 9.2

# /*------------------------------------------------------------------*/
# /*    Based on the SAS macro program, POWERLOG.SAS (version 1.1)    */
# /*    Created 1998, revised 2002                                    */
# /*    SAS author:  Michael Friendly   <FRIENDLY@YorkU.CA>           */
# /*------------------------------------------------------------------*/

# /*------------------------------------------------------------------*/
# /*       p1 -- estimated probability of a 1 at the mean             */
# /*             of the predictor                                     */
# /*       p2 -- estimated probability of a 1 at an X-value           */
# /*             equal to the X-mean plus 1 sd                        */
# /*      rsq -- squared multiple correlation of the predictor        */
# /*             with all other predictors.  Use rsq=0 for a          */
# /*             1-predictor model.                                   */
# /*    alpha -- significance level of test of factor A (default=.05) */
# /*------------------------------------------------------------------*/

"""
    powerlog(p1::Float64,p2::Float64; alpha::Float64 = 0.05, rsq::Real = 0, graph = false, help = false)

Perform power analysis for logistic regression designs with a
continuous predictor variable.  Multiple predictor models can be accomodated using
the `rsq` option. This is a port of Stata powerlog.ado, which was originally ported
from the SAS macro program, POWERLOG.SAS.

## Options
- `p1` = estimated probability of an event (1) at the mean of the predictor
- `p2` = estimated probability of a `1` at the mean + 1 sd of the predictor (required)
- `rsq` = squared multiple correlation of the predictor with all other predictors. Use rsq = 0 for a 1-predictor model (default: 0.0).
- `alpha` = the alpha level, probability of type I error (default: 0.05)
- `graph` = set to `true` to generate a graph of the power analysis (default: false)
- `help` = set to `true` to display explanation of terms

## Examples

```jldoctest
julia> powerlog(.25,.35,alpha = 0.01, help = true)

Logistic regression power analysis
One-tailed test: alpha = 0.01, p1 = 0.25, p2 = 0.35, rsq = 0, odds ratio = 1.6153846153846154

 power        n
  0.60      192
  0.65      211
  0.70      232
  0.75      256
  0.80      284
  0.85      319
  0.90      365
  0.95      439

Explanation of terms
p1  -- the probability that the response variable equals 1
         when the predictor is at the mean
p2  -- the probability that the response variable equals 1
         when the predictor is one standard deviation above the mean
rsq -- the squared mulitple correlation between the predictor
         variable and all other variables in the model

julia> powerlog(.25,.35,rsq=0.4)

Logistic regression power analysis
One-tailed test: alpha = 0.05, p1 = 0.25, p2 = 0.35, rsq = 0.4, odds ratio = 1.6153846153846154

 power          n
 0.60         173
 0.65         196
 0.70         223
 0.75         253
 0.80         290
 0.85         335
 0.90         397
 0.95         498

```

## Acknowledgements
Based on (copied from) the SAS macro program, POWERLOG.SAS (version 1.1),
Created 1998, revised 2002, by  Michael Friendly  <FRIENDLY@YorkU.CA>,
and STATA program, powerlog.ado, by Philip B. Ender <ender@ucla.edu>,
UCLA Academic Technology Services.

"""
function powerlog(p1::Float64,p2::Float64; alpha::Float64 = 0.05, rsq::Real = 0, graph = false, help = false)

    pd = p2 - p1
    l1 = p1/(1-p1)
    l2 = p2/(1-p2)
    θ = l2 / l1
    or = θ
    λ = log(θ)
    λ2 = λ^2
    za = quantile(Normal(),1-alpha)

    # start output
    println("\nLogistic regression power analysis")
    println("One-tailed test: alpha = ",alpha,", p1 = ",p1,", p2 = ",p2,", rsq = ",rsq,", odds ratio = ",or)
    println("\n power          n")

    δ = (1 + (1 + λ2)*exp(5 * λ2/4))/(1 + exp(-1*λ2/4))

    pwr = zeros(Float64,8)
    nn = zeros(Int64,8)

    i = 1
    for power = 0.6:.05:.95
        zb = quantile(Normal(),power)

        N = ((za + zb*exp(-1 * λ2/4))^2 * (1 + 2*p1*δ))/(p1*λ2)
        N /= (1 - rsq)

        @printf(" %3.2f  %10.0f\n",power,N)

        # save for a graph
        pwr[i] = power
        nn[i] = ceil(Int64,N)
        i += 1
    end

    if help == true
        println("\nExplanation of terms")
        println("p1  -- the probability that the response variable equals 1")
        println("         when the predictor is at the mean")
        println("p2  -- the probability that the response variable equals 1")
        println("         when the predictor is one standard deviation above the mean")
        println("rsq -- the squared mulitple correlation between the predictor")
        println("         variable and all other variables in the model")
    end

    if graph == true
        gr()
        len = length(nn)
        n_lower = nn[1]
        n_upper = nn[len]
        xlim_lower = nn[len] / 80
        plot(nn,pwr,
            linecolor = :blue,
            marker=(2,.5,:circle,:blue),
            title = "Logistic regression power analysis",
            yaxis = ("Power",(-0.03,1.03),0.0:0.2:1.0),
            xaxis = ("Sample Size",(n_lower - xlim_lower,n_upper + xlim_lower)),
            legend=false)
        vline!([nn[5]],line=(1,:dot,.8,:red))
    end
end


end # module
