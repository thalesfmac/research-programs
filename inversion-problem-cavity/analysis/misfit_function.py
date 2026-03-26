"""Módulo com funções usadas nas análises"""

from typing import Literal

import numpy as np
from numpy.typing import NDArray
from scipy.integrate import simpson, trapezoid


def integrate(
    x: NDArray[np.float64],
    y: NDArray[np.float64],
    axis: int = -1,
    method: Literal["trapezoid", "simpson", "sum"] = "trapezoid",
):
    """
    Integra um conjunto de dados `y(x)` ao longo de um eixo especificado,
    usando diferentes métodos numéricos.

    Parâmetros
    ----------
    x : NDArray[np.float64]
        Array unidimensional com os pontos da variável independente.
        Deve estar ordenado e ter o mesmo tamanho que `y` ao longo do
        eixo de integração.
    y : NDArray[np.float64]
        Array com os valores da função a serem integrados. Pode ser
        multidimensional, mas deve ser compatível com `x` ao longo do
        eixo dado por `axis`.
    axis : int, opcional
        Eixo ao longo do qual a integração será realizada. O padrão é -1
        (último eixo).
    method : {"trapezoid", "simpson", "sum"}, opcional
        Método de integração numérica:

        - "trapezoid": regra do trapézio (`scipy.integrate.trapezoid`);
        - "simpson": regra de Simpson (`scipy.integrate.simpson`);
        - "sum": soma discreta `sum(y) * Δx`, assumindo espaçamento
          uniforme em `x`.

    Retorna
    -------
    integral : np.ndarray ou np.float64
        Valor da integral ao longo do eixo especificado. O tipo e o
        formato seguem o retorno das rotinas do SciPy ou da soma, a
        depender do método escolhido.

    Levanta
    -------
    ValueError
        Se `method` não for um dos valores suportados
        {"trapezoid", "simpson", "sum"}.

    Notas
    -----
    - Para `method="sum"`, assume-se que o espaçamento em `x` é uniforme,
      e `Δx` é estimado como `np.diff(x)[0]`.
    - Para `method="simpson"`, aplicam-se as mesmas restrições da rotina
      `scipy.integrate.simpson` (por exemplo, número de pontos pode
      precisar ser ímpar, dependendo do caso).
    """
    if method == "trapezoid":
        integral = trapezoid(x=x, y=y, axis=axis)
    elif method == "simpson":
        integral = simpson(x=x, y=y, axis=axis)
    elif method == "sum":
        delta_e = np.diff(x)[0]  # assume espaçamento uniforme
        integral = np.sum(y, axis=axis) * delta_e
    else:
        raise ValueError(f"Unknown integration method: {method}")

    return integral


def misfit_function(
    energies: NDArray[np.float64],
    input_gamma: NDArray[np.float64],
    ca_gamma: NDArray[np.float64],
    e_menos: float,
    e_mais: float,
    energies_axis: int = -1,
    method: Literal["trapezoid", "simpson", "sum"] = "trapezoid",
):
    if e_menos >= e_mais:
        raise ValueError(f"{e_mais=} has to be greater than {e_menos=}")
    mask = (energies >= e_menos) & (energies <= e_mais)

    integrand = (input_gamma - ca_gamma) ** 2

    # Aplica a máscara ao último eixo (energias)
    masked_integrand = np.compress(mask, integrand, axis=energies_axis)
    masked_energies = np.compress(mask, energies)

    integral = integrate(
        x=masked_energies,
        y=masked_integrand,
        axis=energies_axis,
        method=method,
    )

    integral /= e_mais - e_menos
    return integral  # type: ignore


def deltai_function(
    energies: NDArray[np.float64],
    ca_gamma: NDArray[np.float64],
    ca_gamma_ref: NDArray[np.float64],
    e_menos: float,
    e_mais: float,
    energies_axis: int = -1,
    method: Literal["trapezoid", "simpson", "sum"] = "trapezoid",
):
    if e_menos >= e_mais:
        raise ValueError(f"{e_mais=} has to be greater than {e_menos=}")
    mask = (energies >= e_menos) & (energies <= e_mais)

    integrand_0 = ca_gamma
    integrand_1 = np.abs(ca_gamma_ref - ca_gamma)
    integrand_2 = ca_gamma_ref

    # Aplica a máscara ao último eixo (energias)
    masked_energies = np.compress(mask, energies)
    masked_integrand_0 = np.compress(mask, integrand_0, axis=energies_axis)
    masked_integrand_1 = np.compress(mask, integrand_1, axis=energies_axis)
    masked_integrand_2 = np.compress(mask, integrand_2, axis=energies_axis)

    integral_0 = integrate(
        x=masked_energies,
        y=masked_integrand_0,
        axis=energies_axis,
        method=method,
    )
    integral_1 = integrate(
        x=masked_energies,
        y=masked_integrand_1,
        axis=energies_axis,
        method=method,
    )
    integral_2 = integrate(
        x=masked_energies,
        y=masked_integrand_2,
        axis=energies_axis,
        method=method,
    )

    delta = integral_1 / np.abs(integral_2)

    print("Integral ca: ", integral_0)
    print("Integral dif:", integral_1)
    print("Integral ref:", integral_2)
    # return integral_1  # type: ignore
    return delta  # type: ignore
