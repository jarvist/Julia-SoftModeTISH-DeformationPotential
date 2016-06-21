# MAPI-demo.jl
# Demonostration / scratchpad, orientated towards MAPI zone boundary acoustic tilts

using PyPlot
pygui(true)

push!(LOAD_PATH,"./") # Temporary versions of modules in PWD
using SoftModeTISH 

fig=figure(figsize=(12,8)) # Larger than normal figures, for nicer full screen graphs for talks

#TISH(r->r^2) # Harmonic well
#TISH(r->0) # Infinite well potential.
#TISH(r->3*r^2+r)
#TISH(r->-0.1/abs(r)) # Atom like - 1/r Coulomb potential
#TISH(r->20*r^2-20*abs(r),100,6) # Pretty crummy mexican hat
TISH(r->10E-2*r^4-8E-2*r^2,99,5) # Double well potential; for a 'soft mode' quantum harmonic oscillator phonon instability
#TISH(r->10E-2*r^4-8E-2*r^2,99,0) # Just visualise the above potential (for figure for talk)
#TISH(r->10E-2*r^4,99,20) # Single well potential


N=99
evals,evecs=TISH(r->-1E-3*(abs(r+0.5))^-1-1E-3*(abs(r-0.5))^-1,N,4) # Two minima, sort of like a hydrogen molecule

totaldensity=0.0
@printf("Bose Einstein weight matrix: ")

fig=figure() # new figure please
xlabel("Q")
ylabel("Bose-Einstein weighted contribution density")

for i in 1:length(evals)
#    println(BE(evals[i],evals[1]) )
    T=300
    BEweight=BE(evals[i],evals[1]-kBeV*T,T)
    # alpha set to kbT below the lowest energy level, ~ unitary summation
    @printf("eval=%f eV BE=%f ",evals[i],BEweight)

    # Plot of weighted ψ^2 probability densities 
    plot(BEweight * evecs[:,i].^2)
    # Sum up density
    totaldensity+= (BEweight * evecs[:,i].^2) # Density of this structure
end

fig=figure() # new figure please
plot(totaldensity) # curve of total density
fill_between(0:N,0,totaldensity,color="blue",alpha=0.3) # nice filled curve, partially transparent
xlabel("Q")
ylabel("Probability Density Function")

sum(totaldensity)

#OK; now we want to calculate some kind of electron phonon coupling, so let's
#try convolving the deformation potential (tabulated) with the density from
#these thingies

N=99
evals,evecs=TISH(r->10E-2*r^4-8E-2*r^2,N,5) # Double well potential; for a 'soft mode' quantum harmonic oscillator phonon instability

fig=figure() # new figure please
xlabel("Q")
ylabel("Probability Density Function")


WavefunctionDensity=evecs[:,1].^2
plot(WavefunctionDensity,color="magenta")    # Potential energy curve

DeformationPotential(Q)=16E-2*Q^2 # Approximate fit to Lucy's data; 160 meV at the extrema, quadratic
DeformationTabulated=[DeformationPotential(r)::Float64 for r in -1.0:2/N:1.0]
plot(DeformationTabulated,color="black")    # Q-resolved deformation potential

EPhCouple=WavefunctionDensity.*DeformationTabulated
plot(EPhCouple,color="green")    # Q-resolved coupling

@printf("E-Ph coupling: %f eV\n",sum(EPhCouple))

### All singing

N=500
dx=1E2/(N-1)
evals,evecs=TISH(r->20E-2*r^4-10E-2*r^2,N,10,dx) # Double well potential; for a 'soft mode' quantum harmonic oscillator phonon instability
#evals,evecs=TISH(r->10E-2*r^4-8E-2*r^2,N,5)
#evals,evecs=TISH(r->10E-2*r^4,N,20) # Single well potential

fig=figure() # new figure please
xlabel("Q")
ylabel("E-Ph coupling")
title("Temperature dependent E-ph coupling")

# Iterate over these temperatures
for T in collect(1:50:1000)  #[1,50,100,150,200,300,600,10000] #collect(1:300) 
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
