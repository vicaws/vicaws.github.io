---
layout: distill
title: Edge waves
date: 2017-11-04 12:00:00-0000
description:
tags: maths, physics
categories: maths
giscus_comments: true

authors:
  - name: Victor Wang
    affiliations:
      name: InFoMM, Oxford

toc:
  - name: Stokes waves model
  - name: Linearisation
  - name: Solution and discussion
---

The waves on the sea surface at the beach can be modelled as the Stokes waves propagating on a sloping bottom. In the following sections, we have written down the governing equations for the surface waves, and then derive the dispersion relations for these waves.


## Stokes waves model

The waves on the sea are typical stoke waves, the waves at the interface between a inviscid, irrotational incompressible fluid and a gas, for which the modelling assumes constant pressure and no shear stress on the interface, and that the only external force on the system is gravity acting on the fluid.

Take the beach to be at $x=0$ abd the sea to extend in $0<x \rightarrow \infty$ with $-\infty < y < \infty$. The undisturbed surface of the sea is $z = 0$ and the sea bed slopes away at an angle $\beta$. Consider an inviscid, irrotational incompressible fluid occupying the region $-x \tan \beta < z < h(x,y,t)$, where the bottom surface ($z = -x \tan \beta$) is rigid and impermeable, and the top surface ($z = h(x,y,t)$) is a free boundary. We describe the properties of the fluid by equation (\ref{eq:Irrotationality}) and (\ref{eq:Incompressibility}).

\begin{align}
    \text{Irrotationality} \quad & \exists \phi \text{ s.t. } \mathbf{u} = \nabla \phi,
    \label{eq:Irrotationality}
\end{align}
\begin{align}
    \text{Incompressibility} \quad & \nabla \cdot \mathbf{u} = 0,
    \label{eq:Incompressibility}
\end{align}
where the fluid velocity $\mathbf{u}$ can be expressed as the gradient of a scalar function velocity potential $\phi(x,y,z,t)$.

Now consider boundary conditions. The sea bed is impermeable, so that

\begin{equation}
\mathbf{u} \cdot \mathbf{n} = 0, \quad \text{on } z = -x \tan \beta
\end{equation}
where the unit normal vector is $\mathbf{n} = (\sin \beta, 0, \cos \beta)$. On the free surface, we have the kinetic condition and the dynamic condition: on $ z = h(x,y,t)$,

$$
\text{Kinetic:} \quad \Big( \frac{\partial}{\partial t} + \mathbf{u} \cdot \nabla \Big) (z-h) = 0,
$$

\begin{equation}
\text{Dynamic:} \quad \frac{\partial \phi}{\partial t} + \frac{1}{2} |\nabla \phi|^2 = -gh,
\label{eq:SurfBond}
\end{equation}

## Linearisation

We replace $\mathbf{u}$ by $\nabla \phi$ in all equations, and linearise them by considering small perturbations to a flat free surface and a static fluid:

\begin{equation}
\phi = \delta \bar{\phi}, \quad h = \delta \bar{h}.
\end{equation}

Keeping only leading order terms in $\delta$, and dropping the overbar notations, we have

\begin{equation}
    \nabla^2 \phi = 0, \quad \text{in } -x \tan \beta < z < 0,
    \label{eq:LinearGeneral}
\end{equation}
\begin{equation}
    \tan \beta \frac{\partial \phi}{\partial x} + \frac{\partial \phi}{\partial z} = 0, \quad \text{on } z = -x \tan \beta,
    \label{eq:LinearBottom}
\end{equation}
\begin{equation}
    -\frac{\partial h}{\partial t} + \frac{\partial \phi}{\partial z} = 0, \quad \text{on } z = 0,
    \label{eq:LinearSurfKinetic}
\end{equation}
\begin{equation}
    \frac{\partial \phi}{\partial t} + gh = 0, \quad \text{on } z = 0.
    \label{eq:LinearSurfDynamic}
\end{equation}

## Solution and discussion

Too find solution from equations (\ref{eq:LinearGeneral})--(\ref{eq:LinearSurfDynamic}), we assume the waves have form $\phi(x,y,z,t) = f(y-ct) F(x,z)$. Then equation (\ref{eq:LinearGeneral}) gives

\begin{equation}
\label{eq:SepVariables}
\frac{\partial_{xx} F + \partial_{zz} F}{F} = - \frac{\partial_{yy} f}{f} = \lambda^2,
\end{equation}

where $\lambda$ is an arbitrary constant. The equation for $\phi$ shows wave phenomena, so that we attempt for a wave solution

\begin{equation}
f(y-ct) = C_1 e^{i\lambda (y-ct)},
\end{equation}

with a constant $C_1$. Also we seek for a decayed wave solution for the height of the sea with the following form

\begin{equation}
h(x,y,t) = C_2 e^{i \lambda (y-ct)} e^{-ax},
\end{equation}

in which the last term means the waves decay when $x \rightarrow \infty$ with a position constant $a$.

We further assume, by separation of variables, $F(x,z) = G(x)H(z)$. Then equation (\ref{eq:LinearSurfDynamic}) gives

\begin{equation}
G(x) \cdot H(z) \cdot f \cdot (-i \lambda c) = -g C_3 \cdot f e^{-ax}.
\end{equation}

Picking $x$-related terms on both sides of the equation, we find

\begin{equation}
G(x) = C_4 e^{-ax},
\end{equation}
for some constant $C_4$.

Now equation (\ref{eq:SepVariables}) can be written into

\begin{equation}
\frac{H''}{H} = \lambda^2 - a^2.
\end{equation}

We can verify that $H$ is not complex by checking the boundary condition on the sea surface. In fact, equation (\ref{eq:LinearSurfKinetic}) and (\ref{eq:LinearSurfDynamic}) give

\begin{equation}
\left. \frac{H'}{H} \right|_{z=0}= \frac{\lambda^2 c^2}{g}
\label{eq:FinalSurf}
\end{equation}

Hence the solution of $H$ is in the form of

\begin{equation}
H(z) = C_5 e^{kz} + C_6 e^{-kz},
\end{equation}

where $k = \sqrt{\lambda^2 - a^2}$, and $C_5$, $C_6$ are constants. Using the boundary condition on the sea bed, namely equation (\ref{eq:LinearBottom}), we have

\begin{equation}
\label{eq:FinalSeaBed}
-a \tan \beta + k \frac{C_5 e^{kz} - C_6 e^{-kz}}{C_5 e^{kz} + C_6 e^{-kz}} = 0, \text{ at } z = -x \tan \beta.
\end{equation}
Equation (\ref{eq:FinalSeaBed}) holds true only when $k = a \tan \beta$ and $C_6 = 0$. Along with the definition of $k$, we obtain

\begin{equation}
k = \lambda \sin \beta,
\end{equation}

and therefore

\begin{equation}
\label{eq:HDerivative}
\frac{H'}{H} = \lambda \sin \beta.
\end{equation}

Note that the wave frequency $\omega$ is denoted as $\lambda c$ in all the abovementioned equations. Equation (\ref{eq:FinalSurf}) and (\ref{eq:HDerivative}) lead to the dispersion relation between wave frequency $\omega$ and wave number $\lambda$.

\begin{equation}
\label{eq:DispersionRelation}
\omega^2 = g \lambda \sin \beta.
\end{equation}

Recall that the dispersion relation for the classical Stokes waves is given by

\begin{equation}
\omega^2 = g \lambda \tanh (\lambda H),
\end{equation}

where $H$ is the constant mean depth of the sea. Consider a sea of infinite depth, i.e. $H \rightarrow \infty$, the classical governing equation becomes

\begin{equation}
\omega^2 = g \lambda,
\end{equation}

which agrees with the relation (\ref{eq:DispersionRelation}) when $\beta = \frac{\pi}{2}$.
