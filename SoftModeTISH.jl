module SoftModeTISH

using PyPlot 

export TISH,TISHplot, BE, kBeV

# Following: http://www.cond-mat.de/teaching/DFT/qm1d.html
# # 1D numeric Shrodinger equation solver, by discretisation to an Eigenvalue problem
# #  $f''(x_i) \approx (f(x_{i-1})-2*f(x_i)+f(x_{i+1}))/dx^2$

# V -> anonymous function generating the potential
# N -> number of points in 
# dx -> Kinetic Energy operator; goes on the offdiagonal elements in real-space
function TISH(V,N=99,n=3,dx=1E2/(N-1))
   
    # PE terms on the trace
    diagonal = [(2.0/dx^2 + V(r))::Float64 for r in -1.0:2/N:1.0]
    
    # KE terms on the tridiagonals
    updiagonal = [(-1/dx^2)::Float64 for r in 1:N]
    H =diagm(diagonal,0) + diagm(updiagonal,1) + diagm(updiagonal,-1)

    # And solve with dense eigensolvers
    evals,evecs=eig(H)
   
    return evals,evecs
end

# n -> number of states to plot (all N states are calculated)
function TISHplot(V,evals,evecs, n=3)
    N=length(evals)
    xlabel("Q")
    ylabel("Energy (eV) (+ Psi^2)")
    plot([V(r) for r in -1.0:2/N:1.0],color="black")    # Potential energy curve
    
    # This many eigenenergies
    for i in 1:n
        # Ψ ; the wavefunction, offset by the eigenvalue
        plot(1E-2.*evecs[:,i]+evals[i])
        
        # Ψ^2 , the Prob. density, plotted grey, offset by the eigenvalues
        plot(1E-1.*evecs[:,i].^2+evals[i],color="grey")
        # Ψ^2, the prob density, filled curve in semi-tranparent grey, offset + to the eigenvalues
        fill_between(1:N,evals[i],evals[i]+1E-1.*evecs[:,i].^2,color="grey",alpha=0.3)
    end
end 

# Right - very well! We have a 1D TISH solver.  

# Now we're going to set up a function to generate a Bose-Einstein
# distribution, to be able to calculate the occupation of these states in a
# Quantum Harmonic Oscillator (i.e. phonon
# populations)

kBeV=8.6173324E-5 # in units of eV

#T=300           # Temperature; Kelvin
#μ=0.0           # Chemical potential; eV
function BE(E,μ=0.0,T=300.0)
    # System specifics + derived quantities
    β=1/(kBeV*T)      # Thermodynamic Beta; units
    1/(exp((E-μ)*β)-1)
end

# So now, let's get our `evals` and `evecs` from the TISH solver, then
# integrate over the set of eigenvalues, multiplying each probability density
# function ($\Psi^2$) with a Bose-Einstein (BE) occupation factor.
#
# We'll set the internal $\mu$ in the BE distribution to a value kbT below the
# lowest energy level to make it much less infinity, and roughly integrate out
# to give ~1, but then explicitly normalise the total density to exactly 1. I'm
# not quite sure what value should be used for $\mu$, it acts as
# a normalisation constant for the DoS. I assume it is explicitly a chemical
# potential, but as the population of Phonons changes, this confuses me.

end
