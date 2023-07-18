tmp<fv::convectionScheme<scalar>> mvConvection
(
    fv::convectionScheme<scalar>::New
    (
        mesh,
        fields,
        rhoPhi,
        mesh.divScheme("div(phi,Yi_h)")
    )
);

{
    volScalarField Yt(0.0*Y[0]);
        reaction->correct();
        Qdot = reaction->Qdot();
	    forAll(Y, i)
	    {
		if (i != inertIndex && composition.active(i))
		{
		    volScalarField& Yi = Y[i];
			
		    fvScalarMatrix YiEqn
		    (
		        fvm::ddt(rho, Yi)
		      + mvConvection->fvmDiv(rhoPhi, Yi)
		      + transport.Yflux(Yi)
		     ==
		        reaction->R(Yi)
		    );

		    YiEqn.relax();
		    YiEqn.solve(mesh.solver("Yi"));

		    Yi.max(0.0);
                    Yi.correctBoundaryConditions();
		    Yt += Yi;
		}
	    }
   

    Y[inertIndex] = scalar(1) - Yt;
    Y[inertIndex].max(0.0);
}
