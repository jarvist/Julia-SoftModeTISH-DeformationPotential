# MAPI-demo.jl
# Demonostration / scratchpad, orientated towards MAPI zone boundary acoustic tilts

using PyPlot
pygui(true)

push!(LOAD_PATH,"./") # Temporary versions of modules in PWD
using SoftModeTISH 

### Global variables 
N=1500
QMax=150
#QMax=60 # to compare to Lucy's fit


# Function to wrap everything that was in a script
function publicationimages(modename,V)
#dx=2E2/(N-1) # best guess from formula
dx=6.11e3/(N-1) # This value for the KE gets agreement with John's Fourier Method
evals,evecs=TISH(V,N,dx,QMax) # calculate with module 
println(evals[1]," ",evals[3]," ",evals[4]) # Compare eigenvalues to John's code

fig=figure()
#fig=figure(figsize=(12,8)) # Larger than normal figures, for nicer full screen graphs for talks

TISHplot(V,evals,evecs,200,QMax) # Plot eigenmodes

PyPlot.ylim((-0.04,0.04))
PyPlot.xlim((-70.0,70.0))

PyPlot.tight_layout()
PyPlot.savefig(string(modename,"-Modes.png"), format="png", dpi=600)

print(evals[1],evals[3],evals[4])
#evals,evecs=TISH(r->10E-2*r^4-8E-2*r^2,N,5)
#evals,evecs=TISH(r->10E-2*r^4,N,20) # Single well potential

fig=figure()
#fig=figure(figsize=(7.5/2.54, 16.0/2.54)) # new figure please
# Standardising against Jonathan's Matplotlib settings
xlabel(L"$Q_0$ [amu$^{\frac{1}{2}}$ $\AA$]")
#PyPlot.rc('font', **{ 'family' : 'serif', 'size' : '8', 'serif' : 'Times New Roman' })
#PyPlot.rc('lines', **{ 'linewidth' : 0.5 })

#xlabel("Q")
ylabel("DoS")
#title("Nuclear Coordinate Density of State")


DeformationPotential(Q)=160E-3*(Q/70.0)^2 # Vaguely fitted to Lucy's data; quadratic form
#    DeformationPotential(Q)=16E-2*abs(Q) # Linear form
#    DeformationPotential(Q)=1.0 # constant
DeformationTabulated=[DeformationPotential(Q)::Float64 for Q in -QMax:(2.0*QMax)/N:QMax]
#    plot(DeformationTabulated,color="pink")    # Q-resolved deformation potential


# Iterate over temperature and generate PDF figures
for T in [50,150,300,3000] #collect(1:10:1000)  #[1,50,100,150,200,300,600,10000] #collect(1:300) 
    totaldensity=BEWeightedDensity(evals,evecs,T)
    plot(collect(-QMax:2*QMax/N:QMax),totaldensity,label=@sprintf("%d K",T))    # Density profiles 
    fill_between(collect(-QMax:2*QMax/N:QMax),0,totaldensity,color="grey",alpha=0.15) # Add a filled curve in grey for PDF, with light alpha

    EPhCouple=totaldensity.*DeformationTabulated
#    plot(EPhCouple,label=@sprintf("%d K",T))    # Q-resolved coupling
    
#    fig=figure() # new figure please
#    plt[:hist](EPhCouple,100) 
 
    @printf("T: %d E-Ph: %f eV\n",T,sum(EPhCouple))
    
    # Plot this together in a vaguely pleasant way
#    plot(totaldensity,label=@sprintf("%d K",T)) # Plot the PDF, with a label
end
legend(loc="upper right",fancybox="true") # Create a legend of all the existing plots using their labels as names

PyPlot.tight_layout()
PyPlot.savefig(string(modename,"-PDF.png"), format="png", dpi=300)

for T in collect(1:10:1000)  
    totaldensity=BEWeightedDensity(evals,evecs,T)
    EPhCouple=totaldensity.*DeformationTabulated
#    plot(EPhCouple,label=@sprintf("%d K",T))    # Q-resolved coupling
    
    @printf("T: %d E-Ph: %f eV\n",T,sum(EPhCouple))
end
end
# End of function wrapping


# General Mexican-hat double-well potential
#V=r->20E-2*r^4-10E-2*r^2

#R1 mode
modename="R1"
V=Q-> -3.96276195741e-05*Q^2 + 1.05036527538e-08*Q^4 + 1.50096872115e-15*Q^6 + 7.23239952095e-19*Q^8
publicationimages(modename,V) #calculates, generates figures

#M1 mode
modename="M1"
V=Q-> 8.15357577448e-18*Q + -2.6032640735e-05*Q^2 + 1.44432132845e-21*Q^3 + 8.98657195055e-09*Q^4 + -3.18891361563e-26*Q^5 + -5.16891155313e-14*Q^6 + 2.72146169582e-30*Q^7 + 6.92418668864e-19*Q^8
publicationimages(modename,V)
