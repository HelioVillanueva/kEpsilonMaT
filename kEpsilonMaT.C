/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2021 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "kEpsilonMaT.H"
#include "fvModels.H"
#include "fvConstraints.H"
#include "bound.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class BasicMomentumTransportModel>
tmp<volScalarField::Internal> kEpsilonMaT<BasicMomentumTransportModel>::Matk
(
    const volScalarField& Cp,
    const volScalarField& Cv,
    const volScalarField& p
) const
{
    // specific heat ratio
    const volScalarField::Internal gamma = Cp/Cv;
    // acoustic velocity m/s
    const volScalarField::Internal a = sqrt(gamma*p/this->rho_);

    return sqrt(2*k_.v())/a;
}


template<class BasicMomentumTransportModel>
tmp<volScalarField> kEpsilonMaT<BasicMomentumTransportModel>::fMat
(
    const volScalarField& Mat
) const
{
    const volScalarField Matm1 = Mat - Mat0_;
    volScalarField hs = Matm1;

    forAll(hs,i)
    {
        hs[i] = max(Matm1[i], 0.0)/Matm1[i];
    }

    return ((sqr(Mat) - sqr(Mat0_))*hs);
}

template<class BasicMomentumTransportModel>
tmp<volScalarField> kEpsilonMaT<BasicMomentumTransportModel>::rCmu
(
    const volScalarField& T,
    const volScalarField& Cp,
    const volScalarField& Cv,
    const volScalarField& p
)
{
    // specific heat ratio
    const volScalarField gamma = Cp/Cv;
    // acoustic velocity m/s
    const volScalarField a = sqrt(gamma*p/this->rho_);
    // Mach number
    const volScalarField Ma = mag(this->U_)/a;
    // Turbulent Mach number
    const volScalarField Mat = sqrt(2*k_)/a;
    // total temperature   
    const volScalarField Ttot = T*(1.0 + (gamma - 1)*sqr(Ma)/2);

    const volScalarField Tg
    (
        mag(fvc::grad(Ttot))*pow(k_,1.5)/epsilon_/Ttot
    );
    
    const volScalarField CT(1.0 + 1.2*pow(Tg,0.6)/(1.0 + fMat(Mat)));

    return Cmu_/CT;
}


template<class BasicMomentumTransportModel>
void kEpsilonMaT<BasicMomentumTransportModel>::correctNut
(
    const volScalarField& T,
    const volScalarField& Cp,
    const volScalarField& Cv,
    const volScalarField& p
)
{
    this->nut_ = rCmu(T,Cp,Cv,p)*sqr(k_)/epsilon_;
    this->nut_.correctBoundaryConditions();
    fvConstraints::New(this->mesh_).constrain(this->nut_);
}


template<class BasicMomentumTransportModel>
void kEpsilonMaT<BasicMomentumTransportModel>::correctNut()
{
    const volScalarField& T =
        this->mesh_.objectRegistry::template 
        lookupObject<volScalarField>(TName_);
    const volScalarField& Cp =
        this->mesh_.objectRegistry::template 
        lookupObject<volScalarField>(CpName_);
    const volScalarField& Cv =
        this->mesh_.objectRegistry::template
        lookupObject<volScalarField>(CvName_);
    const volScalarField& p =
        this->mesh_.objectRegistry::template
        lookupObject<volScalarField>("p");

    correctNut(T,Cp,Cv,p);
}


template<class BasicMomentumTransportModel>
tmp<fvScalarMatrix> kEpsilonMaT<BasicMomentumTransportModel>::kSource() const
{
    return tmp<fvScalarMatrix>
    (
        new fvScalarMatrix
        (
            k_,
            dimVolume*this->rho_.dimensions()*k_.dimensions()
            /dimTime
        )
    );
}


template<class BasicMomentumTransportModel>
tmp<fvScalarMatrix> kEpsilonMaT<BasicMomentumTransportModel>::epsilonSource() const
{
    return tmp<fvScalarMatrix>
    (
        new fvScalarMatrix
        (
            epsilon_,
            dimVolume*this->rho_.dimensions()*epsilon_.dimensions()
            /dimTime
        )
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicMomentumTransportModel>
kEpsilonMaT<BasicMomentumTransportModel>::kEpsilonMaT
(
    const alphaField& alpha,
    const rhoField& rho,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const viscosity& viscosity,
    const word& type
)
:
    eddyViscosity<RASModel<BasicMomentumTransportModel>>
    (
        type,
        alpha,
        rho,
        U,
        alphaRhoPhi,
        phi,
        viscosity
    ),
    Mat0_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Mat0",
            this->coeffDict_,
            0.1
        )
    ),
    Cmu_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cmu",
            this->coeffDict_,
            0.09
        )
    ),
    C1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "C1",
            this->coeffDict_,
            1.44
        )
    ),
    C2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "C2",
            this->coeffDict_,
            1.92
        )
    ),
    C3_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "C3",
            this->coeffDict_,
            0
        )
    ),
    sigmak_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "sigmak",
            this->coeffDict_,
            1.0
        )
    ),
    sigmaEps_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "sigmaEps",
            this->coeffDict_,
            1.3
        )
    ),
    TName_(IOobject::groupName("T", alphaRhoPhi.group())),
    CpName_(IOobject::groupName("Cp", alphaRhoPhi.group())),
    CvName_(IOobject::groupName("Cv", alphaRhoPhi.group())),

    k_
    (
        IOobject
        (
            IOobject::groupName("k", alphaRhoPhi.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    ),
    epsilon_
    (
        IOobject
        (
            IOobject::groupName("epsilon", alphaRhoPhi.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    )
{
    bound(k_, this->kMin_);
    bound(epsilon_, this->epsilonMin_);

    if (type == typeName)
    {
        this->printCoeffs(type);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicMomentumTransportModel>
bool kEpsilonMaT<BasicMomentumTransportModel>::read()
{
    if (eddyViscosity<RASModel<BasicMomentumTransportModel>>::read())
    {
        Mat0_.readIfPresent(this->coeffDict());
        Cmu_.readIfPresent(this->coeffDict());
        C1_.readIfPresent(this->coeffDict());
        C2_.readIfPresent(this->coeffDict());
        C3_.readIfPresent(this->coeffDict());
        sigmak_.readIfPresent(this->coeffDict());
        sigmaEps_.readIfPresent(this->coeffDict());

        return true;
    }
    else
    {
        return false;
    }
}


template<class BasicMomentumTransportModel>
void kEpsilonMaT<BasicMomentumTransportModel>::correct()
{
    if (!this->turbulence_)
    {
        return;
    }

    // Local references
    const alphaField& alpha = this->alpha_;
    const rhoField& rho = this->rho_;
    const surfaceScalarField& alphaRhoPhi = this->alphaRhoPhi_;
    const volVectorField& U = this->U_;
    volScalarField& nut = this->nut_;
    const Foam::fvModels& fvModels(Foam::fvModels::New(this->mesh_));
    const Foam::fvConstraints& fvConstraints
    (
        Foam::fvConstraints::New(this->mesh_)
    );

    eddyViscosity<RASModel<BasicMomentumTransportModel>>::correct();

    volScalarField::Internal divU
    (
        fvc::div(fvc::absolute(this->phi(), U))()
    );

    tmp<volTensorField> tgradU = fvc::grad(U);
    volScalarField::Internal G
    (
        this->GName(),
        nut()*(dev(twoSymm(tgradU().v())) && tgradU().v())
    );
    tgradU.clear();

    // Update epsilon and G at the wall
    epsilon_.boundaryFieldRef().updateCoeffs();

    // Dissipation equation
    tmp<fvScalarMatrix> epsEqn
    (
        fvm::ddt(alpha, rho, epsilon_)
      + fvm::div(alphaRhoPhi, epsilon_)
      - fvm::laplacian(alpha*rho*DepsilonEff(), epsilon_)
     ==
        C1_*alpha()*rho()*G*epsilon_()/k_()
      - fvm::SuSp(((2.0/3.0)*C1_ - C3_)*alpha()*rho()*divU, epsilon_)
      - fvm::Sp(C2_*alpha()*rho()*epsilon_()/k_(), epsilon_)
      + epsilonSource()
      + fvModels.source(alpha, rho, epsilon_)
    );

    epsEqn.ref().relax();
    fvConstraints.constrain(epsEqn.ref());
    epsEqn.ref().boundaryManipulate(epsilon_.boundaryFieldRef());
    solve(epsEqn);
    fvConstraints.constrain(epsilon_);
    bound(epsilon_, this->epsilonMin_);

    // Turbulent Mach number for k source term
    const volScalarField& Cp =
        this->mesh_.objectRegistry::template 
        lookupObject<volScalarField>(CpName_);
    const volScalarField& Cv =
        this->mesh_.objectRegistry::template 
        lookupObject<volScalarField>(CvName_);
    const volScalarField& p =
        this->mesh_.objectRegistry::template 
        lookupObject<volScalarField>("p");

    const volScalarField::Internal MaT(this->Matk(Cp,Cv,p));

    // Turbulent kinetic energy equation
    tmp<fvScalarMatrix> kEqn
    (
        fvm::ddt(alpha, rho, k_)
      + fvm::div(alphaRhoPhi, k_)
      - fvm::laplacian(alpha*rho*DkEff(), k_)
     ==
        alpha()*rho()*G
      - fvm::SuSp((2.0/3.0)*alpha()*rho()*divU, k_)
      - fvm::Sp(alpha()*rho()*epsilon_()/k_(), k_)
      - fvm::Sp(alpha()*rho()*epsilon_()*sqr(MaT)/k_(), k_)
      + kSource()
      + fvModels.source(alpha, rho, k_)
    );

    kEqn.ref().relax();
    fvConstraints.constrain(kEqn.ref());
    solve(kEqn);
    fvConstraints.constrain(k_);
    bound(k_, this->kMin_);

    const volScalarField& T =
        this->mesh_.objectRegistry::template 
        lookupObject<volScalarField>(TName_);

    correctNut(T,Cp,Cv,p);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// ************************************************************************* //
