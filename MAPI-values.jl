# MAPI-demo.jl
# Demonostration / scratchpad, orientated towards MAPI zone boundary acoustic tilts

using PyPlot
pygui(true)

push!(LOAD_PATH,"./") # Temporary versions of modules in PWD
using SoftModeTISH 

fig=figure(figsize=(12,8)) # Larger than normal figures, for nicer full screen graphs for talks

### All singing

N=1000
dx=1E2/(N-1)
evals,evecs=TISH(r->20E-2*r^4-10E-2*r^2,N,10,dx) # Double well potential; for a 'soft mode' quantum harmonic oscillator phonon instability
#evals,evecs=TISH(r->10E-2*r^4-8E-2*r^2,N,5)
#evals,evecs=TISH(r->10E-2*r^4,N,20) # Single well potential

fig=figure() # new figure please
xlabel("Q")
ylabel("E-Ph coupling")
title("Temperature dependent E-ph coupling")

# Iterate over these temperatures
for T in collect(1:10:1000)  #[1,50,100,150,200,300,600,10000] #collect(1:300) 
    totaldensity=0.0
    for i in 1:length(evals)
        BEweight=BE(evals[i],evals[1]-kBeV*T,T)
        # alpha set to kbT below the lowest energy level, ~ unitary summation
#        if (BEweight>0.001)
#            @printf("T: %03d State: %d : %f eV BE=%f \n",T,i,evals[i],BEweight)
#        end
        
        # Plot of weighted ψ^2 probability densities 
        #plot(BEweight * evecs[:,i].^2)
        # Sum up density
        totaldensity+= (BEweight * evecs[:,i].^2) # Density of this structure
    end
    totaldensity/=sum(totaldensity) # renormalise probability density to ∫ dx =1
    
    DeformationPotential(Q)=16E-2*Q^2 # Vaguely fitted to Lucy's data; quadratic form
#    DeformationPotential(Q)=16E-2*abs(Q) # Linear form
#    DeformationPotential(Q)=1.0 # constant
    
    DeformationTabulated=[DeformationPotential(Q)::Float64 for Q in -1.0:2/N:1.0]
#    plot(DeformationTabulated,color="pink")    # Q-resolved deformation potential

    EPhCouple=totaldensity.*DeformationTabulated
    plot(EPhCouple,label=@sprintf("%d K",T))    # Q-resolved coupling
    
#    fig=figure() # new figure please
#    plt[:hist](EPhCouple,100) 
 
    @printf("T: %d E-Ph: %f eV\n",T,sum(EPhCouple))
    
    # Plot this together in a vaguely pleasant way
#    plot(totaldensity,label=@sprintf("%d K",T)) # Plot the PDF, with a label
#    fill_between(0:N,0,totaldensity,color="grey",alpha=0.05) # Add a filled curve in grey for PDF, with light alpha
end
legend(loc="upper right",fancybox="true") # Create a legend of all the existing plots using their labels as names


#show() # Display any figures produced
