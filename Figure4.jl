##--This Julia File Can Be Used to Generate the Files for All Panels in Figure 4. 
##--This is currently set up for SD=0 & CoV=1 (Panel C). Where to change is indicated by CHANGE HERE
##--This creates 2 files, you will need both. Where the files save is determined by work_dir
using DataFrames, Random, CSV, Distributions, Distributed, Plots, RCall

function namedfile1(df,namefile)

    work_dir="//Users//nikas//Desktop//Example-RiskCorrelate-" #This is where the files will save out, change if on a different computer
    name = work_dir*namefile
    CSV.write(name,df)
end

function AbWane(t, N, protoA)
    #Gives Ab distribution 
    N=100000
    A=protoA[1:N]
    #A=2
    Y=A*128*((t+41)/42)^(-1) #Waning
    VEab=Y.^(-0.5) #This is technically reduction in suseptibility (ie 1-VE)
    for k in 1:N
        if VEab[k]>1
            VEab[k] =1
        elseif VEab[k]<0
            VEab[k] =0
        else
        end
    end
    return Y, VEab
end

function censor(x, m, D)
    #Right Censoring For Cox Model
    for n in 1:m
        if x[n]==0
            x[n] = D
        end
    end
    return x
end

function giveHazards(N, αβ, ββ)
    #Gives underlying susceptibilty (also set up to give beta distribution VE but unused)
    proto=R" rgamma(n=100000, shape=1, rate=400)" #CHANGE HERE: for different CoV
    Risk=proto[1:N]
    VEhetONLY=rand(Beta(αβ,ββ),N)
    for i in 1:N
        if VEhetONLY[i]<0
            VEhetONLY[i]=0
        elseif VEhetONLY[i]>1
            VEhetONLY[i]=1
        else
        end
    end
    return  Risk, VEhetONLY
end

sim_number=1 #How many simulations to run
name = collect(1:sim_number)
name2 = collect((1+sim_number):(2*sim_number))
starttime = time()
Random.seed!(383) #Set Random Seed


Threads.@threads for g in 1:sim_number
    D = 365         #Number of days in the season
    N = 100000      #Population size
    vacc_days= 1    #Which days do you vaccinate
    vacc_freq=0.4   #Vaccine frequency
    init_freq=500.0/100000. #Infection Seeders
    intr_day =1     #Day disease is introduced

    #Disease Parameters/Functions
    γ = 5.0 #Time until recovery
    protoA=R"rlnorm(100000, meanlog=0, sdlog=0)" #CHANGE HERE: for different SD (last number)
    attempt1=giveHazards(N,0.1, 0.1) #Give underlying susceptibility


  
    #--Initializing Vectors
        #Individual Level Information
    vaccinated = zeros(Int64,N,3) #Save information about vaccinated (col 1 vac/not, col 2 when, col 3 protection if all/nothing)
    infected = zeros(Int64,N)  # 1 infected, 0 not
    recovered = zeros(Int64,N) # 1 recovered, 0 not
    infected_time = zeros(Int64,N) #when people get infected
    recovered_time = zeros(Int64,N) #when people recover
    ever_infected  = zeros(Int64,N) #if they were ever previously infected
    trialVE = zeros(Float64, N) #Individual level VE
        #Daily Information
    exposures = zeros(Int64,D,2) # saves new exposures per day
    infections = zeros(Int64,D,2) # saves new infections per day
    dailyve = zeros(Float64,D)    #saves Average VE per day
    additional_vac =0             #in case vaccination is spread
   


    for i in 1:D
        zzz=AbWane(i,100000, protoA) #Find everyone's Ab level that day
        trialVE=zzz[2] #This is an individual's protection based on Ab level
        dailyve[i]=mean(trialVE) #daily mean VE
        acquire=rand(N) #Set random chance to acquire level
        if i in vacc_days
            #Figure out who to vaccinate
            unvaccinated = collect(1:N)[(vaccinated[:,1] .==0) .& (ever_infected .== 0)]
            if length(unvaccinated) >= round(Int,vacc_freq/length(vacc_days)*N) #because of the rounding we typically get 39.99% vaccinated
                vaccinate = sample(unvaccinated,round(Int,vacc_freq/length(vacc_days)*N),replace = false)
                vaccinated[vaccinate,1] .= 1#mark as vaccinated
                vaccinated[vaccinate,2] .= i #save day of vaccination
                vaccinated[vaccinate,3] .= 1# protected (1) or not (0)(this changes if all or nothing)
                additional_vac_add =(vacc_freq/length(vacc_days)*N)-round(Int,vacc_freq/length(vacc_days)*N)#if vaccinating over multiple days, may round down too many times so add in a person whan this is >1
                additional_vac=additional_vac+additional_vac_add 
                if additional_vac >1
                    vaccinate = sample(unvaccinated,round(Int,1),replace = false)
                    vaccinated[vaccinate,1] .= 1#mark as vaccinated
                    vaccinated[vaccinate,2] .= i #save day of vaccination
                    vaccinated[vaccinate,3] .= 1
                    additional_vac=0
                end
            end #Ends Going Through People to Be Vaccinated on a Given Day
            println("Vaccinating", i)
        end#Ends Vaccination Loop

        vac_inf = 0
        unv_inf=0
        vac_exp = 0
        unv_exp=0
        inf_exp=0
        rec_exp=0

        for j in 1:N
            if infected[j]==1 #if infected, check to see if recovered 
                if infected_time[j]+γ==i
                    infected[j]=0
                    recovered[j]=1
                    recovered_time[j]=i
                end
            elseif infected[j]==0 && recovered[j]==0 #if not infected, then try to infect
                if vaccinated[j]==0 #if no vaccinated
                    if acquire[j]<attempt1[1][j]
                        infected[j]=1
                        infected_time[j]=i
                        ever_infected[j]=1
                        unv_inf=unv_inf+1
                    else
                        unv_exp=unv_exp+1
                    end
                elseif vaccinated[j]==1 #if vaccinated
                    if acquire[j]<attempt1[1][j]*trialVE[j]
                        infected[j]=1
                        infected_time[j]=i
                        ever_infected[j]=1
                        vac_inf= vac_inf+1
                    else
                        vac_exp =vac_exp+1
                    end
                end
            elseif recovered[j]==1
            end
            
        
        end
        
        exposures[i,1] = vac_exp
        exposures[i,2] = unv_exp 
        infections[i,1]= vac_inf
        infections[i,2]= unv_inf
    end#Ends Daily Loop

    #Censor
    infected_time=censor(infected_time, N, D)
    
    #Save out information
    id=collect(1:N)
    mydf = DataFrame(id=id,vac_status = vaccinated[:,1], ever_infected = ever_infected, infected_time = infected_time, vac_time = vaccinated[:,2],  RR=attempt1[1], RRVE=attempt1[2], Testing=trialVE)
    namedfile1(mydf,"$(g).csv")
    mydf2 = DataFrame(Exp_Vac=exposures[:,1], Exp_Unv=exposures[:,2], Inf_Vac=infections[:,1], Inf_Unv=infections[:,2],DailyVE = dailyve)
    namedfile1(mydf2,"$(g+sim_number).csv") 


    
end 

timer=time()-starttime




