from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Callable

import numpy as np


colors: dict[str, str] = {
    "red": "#EC4C4C",
    "green": "#70D07C",
    "blue": "#4C9CED",
    "yellow": "#FFC857",
}


@dataclass
class GenerationResult:
    s: Any
    q: Any
    w: Any
    s_prime: Any
    q_prime: Any
    cov_ws: Any
    cov_wq: Any
    cov_wq_est_max: Any
    rho_sq: Any
    rho_sq_est: Any
    rho_wq: Any
    rho_ws: Any
    s_hat: Any
    s_hat_est: Any
    E_ws: Any
    E_wq: Any
    E_wq_est_lin: Any
    E_wq_est_phi_bar: Any
    delta_s_bar: Any
    delta_q_bar: Any
    ratio_s: Any
    ratio_q: Any
    QSC: Any
    QSC_ratio: Any


a: float = 2
b: float = 4


def lin_phi_bar(z: np.ndarray) -> np.ndarray:
    return a * z**2 + b


def lin_phi(z: np.ndarray, s_hat: float) -> np.ndarray:
    s_hat_sq = s_hat**2
    sigma_noise = np.sqrt((s_hat_sq**-1 - 1) * np.var(lin_phi_bar(z)))
    noise = np.random.normal(loc=0, scale=sigma_noise, size=z.shape)
    return lin_phi_bar(z) + noise


def lin_phi_inv(q: np.ndarray) -> np.ndarray:
    return (q - b) / a


def noise(k: float, mu: float, sigma: float) -> Callable[[np.ndarray], np.ndarray]:
    return lambda z: (k * z) + mu + np.random.normal(0, sigma, z.shape[0])


def fitness(k: float, mu: float, sigma: float) -> Callable[[np.ndarray], np.ndarray]:
    return lambda z: (z - mu) * k + np.random.normal(0, sigma, z.shape[0])


def logistic_fitness(k: float, s_opt: float) -> Callable[[np.ndarray], np.ndarray]:
    return lambda s: 1 / (1 + np.exp(-k * (s - s_opt)))


def delta(mu: float, sigma: float) -> Callable[[np.ndarray], np.ndarray]:
    return lambda z: np.random.normal(mu, sigma, z.shape[0])


def simulate_generation(
    s: np.ndarray,
    *,
    phi: Callable[[np.ndarray, float], np.ndarray] = lin_phi,
    phi_bar: Callable[[np.ndarray], np.ndarray] = lin_phi_bar,
    s_hat: float = 0.8,
    delta_fn: Callable[[np.ndarray], np.ndarray] = noise(0, 0, 0.1),
    fitness_fn: Callable[[np.ndarray], np.ndarray] = fitness(0.5, 0, 0.2),
) -> GenerationResult:
    w_abs = w = fitness_fn(s)
    w_min = np.min(w)
    if w_min < 0:
        w = w - w_min
    w = w / np.sum(w)

    q = phi(s, s_hat)

    cov_ws = np.cov(w, s, bias=True)[0, 1]
    cov_wq = np.cov(w, q, bias=True)[0, 1]
    rho_sq = np.corrcoef(s, q)[0, 1]
    rho_wq = np.corrcoef(w, q)[0, 1]
    rho_ws = np.corrcoef(w, s)[0, 1]
    s_hat_est = np.sqrt(np.var(phi_bar(s)) / np.var(q))

    delta_s = delta_fn(s)
    E_ws = np.mean(w * delta_s)

    delta_q = phi(s + delta_s, s_hat) - phi(s, s_hat)
    delta_q_phi_bar = phi_bar(s + delta_s) - phi_bar(s)

    E_wq = np.mean(w * delta_q)
    E_wq_est_lin = np.mean(w * (a * delta_s))
    E_wq_est_phi_bar = np.mean(w * delta_q_phi_bar)

    delta_s_bar = cov_ws + E_ws
    delta_q_bar = cov_wq + E_wq

    s_prime = np.random.choice(s, size=s.shape[0], p=w) + delta_s
    q_prime = phi(s_prime, s_hat)

    return GenerationResult(
        s=s,
        q=q,
        w=w_abs,
        s_prime=s_prime,
        q_prime=q_prime,
        cov_ws=cov_ws,
        cov_wq=cov_wq,
        cov_wq_est_max=rho_sq * np.sqrt(np.var(w) * np.var(q)),
        rho_sq=rho_sq,
        rho_sq_est=rho_wq / rho_ws,
        rho_wq=rho_wq,
        rho_ws=rho_ws,
        s_hat=s_hat,
        s_hat_est=s_hat_est,
        E_ws=E_ws,
        E_wq=E_wq,
        E_wq_est_lin=E_wq_est_lin,
        E_wq_est_phi_bar=E_wq_est_phi_bar,
        delta_s_bar=delta_s_bar,
        delta_q_bar=delta_q_bar,
        ratio_s=E_ws / cov_ws if cov_ws != 0 else np.nan,
        ratio_q=E_wq / cov_wq if cov_wq != 0 else np.nan,
        QSC=E_wq / (np.sqrt(np.var(w) * np.var(q))) if np.var(w) * np.var(q) > 0 else np.nan,
        QSC_ratio=(E_wq / np.sqrt(np.var(w) * np.var(q)) / s_hat) if (np.var(w) * np.var(q) > 0 and s_hat != 0) else np.nan,
    )


@dataclass
class GenerationsResult:
    generations: int
    results: Any
    last_result: Any
    s_means: Any
    s_stds: Any
    q_means: Any
    q_stds: Any
    w_means: Any
    w_stds: Any
    rho_sq_means: Any


def simulate_generations(
    *,
    generations: int = 100,
    n: int = 10_000,
    s_hat: float = 0.5,
    **kwargs: Any,
) -> GenerationsResult:
    s = np.random.normal(0, 1, n)

    results: list[GenerationResult] = []
    s_means: list[float] = []
    s_stds: list[float] = []
    q_means: list[float] = []
    q_stds: list[float] = []
    w_means: list[float] = []
    w_stds: list[float] = []

    r: GenerationResult | None = None
    for g in range(generations):
        r = simulate_generation(s_hat=s_hat, s=s, **kwargs)
        s = r.s_prime
        s_means.append(np.mean(r.s_prime))
        s_stds.append(np.std(r.s_prime))
        q_means.append(np.mean(r.q_prime))
        q_stds.append(np.std(r.q_prime))
        w_means.append(np.mean(r.w))
        w_stds.append(np.std(r.w))
        if g % ((generations // 50) + 1) == 0:
            results.append(r)
            print(f"Generation {g}", end="\r")

    assert r is not None
    gr = GenerationsResult(
        generations=generations,
        results=results,
        last_result=r,
        s_means=np.array(s_means),
        s_stds=np.array(s_stds),
        q_means=np.array(q_means),
        q_stds=np.array(q_stds),
        w_means=np.array(w_means),
        w_stds=np.array(w_stds),
        rho_sq_means=np.corrcoef(s_means, q_means)[0, 1] if len(s_means) > 1 else np.nan,
    )
    return gr


