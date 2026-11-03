/* Calculadoras estadísticas · Nicolás Sánchez Álvarez
   Funciones de cálculo puras. Cada calculadora es CALC[slug] = (v, H) => objeto de resultados.
   v: valores de entrada (números, textos de selección, listas o matrices). H: utilidades de formato. */
(function (root) {
  const J = root.jStat || (typeof require !== 'undefined' ? require('jstat') : null);
  const CALC = {};
  const E = msg => { throw new Error(msg); };
  const need = (cond, msg) => { if (!cond) E(msg); };
  const fin = (...xs) => xs.every(Number.isFinite);
  const Phi = z => J.normal.cdf(z, 0, 1);
  const zq = p => J.normal.inv(p, 0, 1);
  const sum = a => a.reduce((s, x) => s + x, 0);
  const mean = a => sum(a) / a.length;
  const vari = a => { const m = mean(a); return sum(a.map(x => (x - m) ** 2)) / (a.length - 1); };
  const lnfact = n => J.gammaln(n + 1);
  const lnchoose = (n, k) => lnfact(n) - lnfact(k) - lnfact(n - k);
  const quantile = (a, p) => { const s = [...a].sort((x, y) => x - y), h = (s.length - 1) * p, lo = Math.floor(h); return s[lo] + (h - lo) * ((s[lo + 1] ?? s[lo]) - s[lo]); };
  const ranks = a => { const idx = a.map((x, i) => [x, i]).sort((p, q) => p[0] - q[0]); const r = new Array(a.length); let i = 0;
    while (i < idx.length) { let j = i; while (j + 1 < idx.length && idx[j + 1][0] === idx[i][0]) j++; const avg = (i + j) / 2 + 1; for (let k = i; k <= j; k++) r[idx[k][1]] = avg; i = j + 1; } return r; };
  const pearson = (x, y) => { const mx = mean(x), my = mean(y); let sxy = 0, sxx = 0, syy = 0; for (let i = 0; i < x.length; i++) { sxy += (x[i] - mx) * (y[i] - my); sxx += (x[i] - mx) ** 2; syy += (y[i] - my) ** 2; } return sxy / Math.sqrt(sxx * syy); };
  const bisect = (f, lo, hi, it = 200) => { let flo = f(lo); for (let i = 0; i < it; i++) { const mid = (lo + hi) / 2, fm = f(mid); if ((fm > 0) === (flo > 0)) { lo = mid; flo = fm; } else hi = mid; } return (lo + hi) / 2; };
  const pTwo = z => 2 * (1 - Phi(Math.abs(z)));
  const pT = (t, df) => 2 * (1 - J.studentt.cdf(Math.abs(t), df));
  const conf = v => { const c = +v.conf; return { c, z: zq(1 - (1 - c) / 2), pc: Math.round(c * 100) }; };

  // Distribuciones no centrales
  const pois = (j, lam) => Math.exp(-lam + j * Math.log(lam) - lnfact(j));
  function ncChiCdf(x, k, lam) {
    if (x <= 0) return 0; if (lam === 0) return J.chisquare.cdf(x, k);
    const h = lam / 2, j0 = Math.floor(h); let s = 0;
    for (let j = j0; j < j0 + 2000; j++) { const w = pois(j, h); s += w * J.chisquare.cdf(x, k + 2 * j); if (w < 1e-14 && j > h) break; }
    for (let j = j0 - 1; j >= 0; j--) { const w = pois(j, h); s += w * J.chisquare.cdf(x, k + 2 * j); if (w < 1e-14) break; }
    return Math.min(1, s);
  }
  function ncFCdf(x, d1, d2, lam) {
    if (x <= 0) return 0; if (lam === 0) return J.centralF.cdf(x, d1, d2);
    const y = d1 * x / (d1 * x + d2), h = lam / 2, j0 = Math.floor(h); let s = 0;
    for (let j = j0; j < j0 + 2000; j++) { const w = pois(j, h); s += w * J.ibeta(y, d1 / 2 + j, d2 / 2); if (w < 1e-14 && j > h) break; }
    for (let j = j0 - 1; j >= 0; j--) { const w = pois(j, h); s += w * J.ibeta(y, d1 / 2 + j, d2 / 2); if (w < 1e-14) break; }
    return Math.min(1, s);
  }
  const nctCdf = (t, df, ncp) => J.noncentralt.cdf(t, df, ncp);
  const ncChiInv = (p, k, lam) => bisect(x => ncChiCdf(x, k, lam) - p, 0, Math.max(100, (k + lam) * 10));

  // Potencia
  const powT = (d, n, alfa, colas, tipo) => {          // tipo: 'ind' (n por grupo) o 'una' (n total)
    const df = tipo === 'ind' ? 2 * n - 2 : n - 1, ncp = tipo === 'ind' ? d * Math.sqrt(n / 2) : d * Math.sqrt(n);
    if (colas === 2) { const tc = J.studentt.inv(1 - alfa / 2, df); return 1 - nctCdf(tc, df, Math.abs(ncp)) + nctCdf(-tc, df, Math.abs(ncp)); }
    const tc = J.studentt.inv(1 - alfa, df); return 1 - nctCdf(tc, df, Math.abs(ncp));
  };
  const powF = (lam, df1, df2, alfa) => 1 - ncFCdf(J.centralF.inv(1 - alfa, df1, df2), df1, df2, lam);
  const powChi = (lam, df, alfa) => 1 - ncChiCdf(J.chisquare.inv(1 - alfa, df), df, lam);
  const nMin = (powf, target, from = 2, max = 1e6) => { let n = from; if (powf(n) >= target) return n;
    let hi = n; while (powf(hi) < target) { n = hi; hi *= 2; if (hi > max) E('El tamaño necesario supera el millón de casos: revisa el efecto esperado.'); }
    let lo = n; while (hi - lo > 1) { const m = Math.floor((lo + hi) / 2); if (powf(m) >= target) hi = m; else lo = m; } return hi; };

  const U = { J, E, need, fin, Phi, zq, sum, mean, vari, lnfact, lnchoose, quantile, ranks, pearson, bisect, pTwo, pT, conf,
    ncChiCdf, ncFCdf, nctCdf, ncChiInv, powT, powF, powChi, nMin };
  root.CALC = CALC; root.CALC_U = U;
  if (typeof module !== 'undefined' && module.exports) module.exports = { CALC, U };
})(typeof window !== 'undefined' ? window : globalThis);

(function (root) {
  const C = root.CALC, U = root.CALC_U, J = U.J, { need, fin } = U;
  const pos = (x, m) => need(x > 0, m);
  const tails = (cdf, x) => ({ izq: cdf, der: 1 - cdf });

  C['distribucion-normal'] = v => { need(fin(v.x, v.mu, v.sigma), 'Rellena todos los campos.'); pos(v.sigma, 'σ debe ser mayor que 0.');
    const F = J.normal.cdf(v.x, v.mu, v.sigma); return { dens: J.normal.pdf(v.x, v.mu, v.sigma), ...tails(F), z: (v.x - v.mu) / v.sigma }; };
  C['normal-intervalo'] = v => { need(fin(v.a, v.b, v.mu, v.sigma), 'Rellena todos los campos.'); pos(v.sigma, 'σ debe ser mayor que 0.'); need(v.a < v.b, 'El límite inferior debe ser menor que el superior.');
    const p = J.normal.cdf(v.b, v.mu, v.sigma) - J.normal.cdf(v.a, v.mu, v.sigma); return { p, fuera: 1 - p }; };
  C['normal-cuantil'] = v => { need(fin(v.p, v.mu, v.sigma), 'Rellena todos los campos.'); need(v.p > 0 && v.p < 1, 'La probabilidad debe estar entre 0 y 1.'); pos(v.sigma, 'σ debe ser mayor que 0.');
    const x = J.normal.inv(v.p, v.mu, v.sigma); return { x, z: (x - v.mu) / v.sigma }; };
  C['z-critico'] = v => { need(v.alfa > 0 && v.alfa < 1, 'α debe estar entre 0 y 1.'); const c = +v.colas;
    const z = U.zq(1 - v.alfa / c); return { z, txt: c === 2 ? `±${z.toFixed(4).replace('.', ',')}` : `${z.toFixed(4).replace('.', ',')}` }; };
  C['area-normal-estandar'] = v => { need(fin(v.z), 'Escribe un valor de z.'); const F = U.Phi(v.z), a = Math.abs(v.z);
    return { izq: F, der: 1 - F, dos: 2 * (1 - U.Phi(a)), entre: 2 * U.Phi(a) - 1, media: U.Phi(a) - 0.5 }; };

  C['distribucion-t'] = v => { need(fin(v.t, v.df), 'Rellena todos los campos.'); pos(v.df, 'Los grados de libertad deben ser mayores que 0.');
    const F = J.studentt.cdf(v.t, v.df); return { dens: J.studentt.pdf(v.t, v.df), ...tails(F), dos: U.pT(v.t, v.df) }; };
  C['t-critico'] = v => { pos(v.df, 'Los grados de libertad deben ser mayores que 0.'); need(v.alfa > 0 && v.alfa < 1, 'α debe estar entre 0 y 1.');
    return { t: J.studentt.inv(1 - v.alfa / +v.colas, v.df) }; };
  C['distribucion-chi'] = v => { need(fin(v.x, v.df), 'Rellena todos los campos.'); pos(v.df, 'Los grados de libertad deben ser mayores que 0.'); need(v.x >= 0, 'χ² no puede ser negativo.');
    const F = J.chisquare.cdf(v.x, v.df); return { dens: J.chisquare.pdf(v.x, v.df), ...tails(F), media: v.df, var: 2 * v.df }; };
  C['chi-critico'] = v => { pos(v.df, 'Los grados de libertad deben ser mayores que 0.'); need(v.alfa > 0 && v.alfa < 1, 'α debe estar entre 0 y 1.');
    return { x: J.chisquare.inv(1 - v.alfa, v.df) }; };
  C['distribucion-f'] = v => { need(fin(v.x, v.df1, v.df2), 'Rellena todos los campos.'); pos(v.df1, 'Los gl deben ser mayores que 0.'); pos(v.df2, 'Los gl deben ser mayores que 0.'); need(v.x >= 0, 'F no puede ser negativo.');
    const F = J.centralF.cdf(v.x, v.df1, v.df2); return { dens: J.centralF.pdf(v.x, v.df1, v.df2), ...tails(F), media: v.df2 > 2 ? v.df2 / (v.df2 - 2) : NaN }; };
  C['f-critico'] = v => { pos(v.df1, 'Los gl deben ser mayores que 0.'); pos(v.df2, 'Los gl deben ser mayores que 0.'); need(v.alfa > 0 && v.alfa < 1, 'α debe estar entre 0 y 1.');
    return { f: J.centralF.inv(1 - v.alfa, v.df1, v.df2) }; };

  const intk = (k, m) => need(Number.isInteger(k) && k >= 0, m);
  C['binomial-probabilidad'] = v => { intk(v.n, 'n debe ser un entero positivo.'); intk(v.k, 'k debe ser un entero no negativo.'); need(v.k <= v.n, 'k no puede superar n.'); need(v.p >= 0 && v.p <= 1, 'p debe estar entre 0 y 1.');
    return { pk: J.binomial.pdf(v.k, v.n, v.p), media: v.n * v.p, var: v.n * v.p * (1 - v.p), de: Math.sqrt(v.n * v.p * (1 - v.p)) }; };
  C['binomial-acumulada'] = v => { intk(v.n, 'n debe ser un entero positivo.'); intk(v.k, 'k debe ser un entero no negativo.'); need(v.k <= v.n, 'k no puede superar n.'); need(v.p >= 0 && v.p <= 1, 'p debe estar entre 0 y 1.');
    const le = J.binomial.cdf(v.k, v.n, v.p), lt = v.k > 0 ? J.binomial.cdf(v.k - 1, v.n, v.p) : 0;
    return { le, lt, ge: 1 - lt, gt: 1 - le }; };
  C['poisson-probabilidad'] = v => { intk(v.k, 'k debe ser un entero no negativo.'); pos(v.lam, 'λ debe ser mayor que 0.');
    return { pk: J.poisson.pdf(v.k, v.lam), media: v.lam, de: Math.sqrt(v.lam) }; };
  C['poisson-acumulada'] = v => { intk(v.k, 'k debe ser un entero no negativo.'); pos(v.lam, 'λ debe ser mayor que 0.');
    const le = J.poisson.cdf(v.k, v.lam), lt = v.k > 0 ? J.poisson.cdf(v.k - 1, v.lam) : 0; return { le, lt, ge: 1 - lt, gt: 1 - le }; };
  C['hipergeometrica'] = v => { [v.N, v.K, v.n, v.k].forEach(x => intk(x, 'Todos los valores deben ser enteros no negativos.'));
    need(v.K <= v.N && v.n <= v.N, 'K y n no pueden superar N.'); need(v.k <= Math.min(v.K, v.n), 'k no puede superar K ni n.');
    const pk = x => (x < Math.max(0, v.n - (v.N - v.K)) || x > Math.min(v.K, v.n)) ? 0 : Math.exp(U.lnchoose(v.K, x) + U.lnchoose(v.N - v.K, v.n - x) - U.lnchoose(v.N, v.n));
    let le = 0; for (let x = 0; x <= v.k; x++) le += pk(x);
    return { pk: pk(v.k), le: Math.min(1, le), gt: Math.max(0, 1 - le), media: v.n * v.K / v.N }; };
  C['geometrica'] = v => { intk(v.k, 'k debe ser un entero mayor o igual que 1.'); need(v.k >= 1, 'k debe ser mayor o igual que 1.'); need(v.p > 0 && v.p <= 1, 'p debe estar entre 0 y 1.');
    return { pk: (1 - v.p) ** (v.k - 1) * v.p, le: 1 - (1 - v.p) ** v.k, media: 1 / v.p }; };
  C['binomial-negativa'] = v => { intk(v.k, 'El número de fracasos debe ser un entero no negativo.'); intk(v.r, 'r debe ser un entero positivo.'); need(v.r >= 1, 'r debe ser al menos 1.'); need(v.p > 0 && v.p <= 1, 'p debe estar entre 0 y 1.');
    const pk = x => Math.exp(U.lnchoose(x + v.r - 1, x) + v.r * Math.log(v.p) + x * Math.log(1 - v.p)); let le = 0; for (let x = 0; x <= v.k; x++) le += pk(x);
    return { pk: pk(v.k), le: Math.min(1, le), media: v.r * (1 - v.p) / v.p }; };
  C['uniforme'] = v => { need(fin(v.a, v.b, v.x), 'Rellena todos los campos.'); need(v.a < v.b, 'a debe ser menor que b.');
    const F = Math.min(1, Math.max(0, (v.x - v.a) / (v.b - v.a))); return { dens: (v.x >= v.a && v.x <= v.b) ? 1 / (v.b - v.a) : 0, izq: F, der: 1 - F, media: (v.a + v.b) / 2, var: (v.b - v.a) ** 2 / 12 }; };
  C['exponencial'] = v => { pos(v.lam, 'λ debe ser mayor que 0.'); need(v.x >= 0, 'x no puede ser negativo.');
    const F = 1 - Math.exp(-v.lam * v.x); return { dens: v.lam * Math.exp(-v.lam * v.x), izq: F, der: 1 - F, media: 1 / v.lam, mediana: Math.log(2) / v.lam }; };
  C['lognormal'] = v => { pos(v.x, 'x debe ser mayor que 0.'); pos(v.sigma, 'σ debe ser mayor que 0.');
    const F = J.lognormal.cdf(v.x, v.mu, v.sigma); return { dens: J.lognormal.pdf(v.x, v.mu, v.sigma), izq: F, der: 1 - F, media: Math.exp(v.mu + v.sigma ** 2 / 2), mediana: Math.exp(v.mu) }; };
  C['distribucion-gamma'] = v => { pos(v.k, 'La forma debe ser mayor que 0.'); pos(v.theta, 'La escala debe ser mayor que 0.'); need(v.x >= 0, 'x no puede ser negativo.');
    const F = J.gamma.cdf(v.x, v.k, v.theta); return { dens: J.gamma.pdf(v.x, v.k, v.theta), izq: F, der: 1 - F, media: v.k * v.theta, var: v.k * v.theta ** 2 }; };
  C['distribucion-beta'] = v => { pos(v.a, 'α debe ser mayor que 0.'); pos(v.b, 'β debe ser mayor que 0.'); need(v.x >= 0 && v.x <= 1, 'x debe estar entre 0 y 1.');
    const F = J.beta.cdf(v.x, v.a, v.b); return { dens: J.beta.pdf(v.x, v.a, v.b), izq: F, der: 1 - F, media: v.a / (v.a + v.b) }; };
  C['t-no-central'] = v => { need(fin(v.t, v.df, v.ncp), 'Rellena todos los campos.'); pos(v.df, 'Los gl deben ser mayores que 0.');
    const F = U.nctCdf(v.t, v.df, v.ncp); return { izq: F, der: 1 - F, dens: J.noncentralt.pdf(v.t, v.df, v.ncp) }; };
  C['f-no-central'] = v => { need(fin(v.x, v.df1, v.df2, v.lam), 'Rellena todos los campos.'); pos(v.df1, 'Los gl deben ser mayores que 0.'); pos(v.df2, 'Los gl deben ser mayores que 0.'); need(v.lam >= 0 && v.x >= 0, 'F y λ no pueden ser negativos.');
    const F = U.ncFCdf(v.x, v.df1, v.df2, v.lam); return { izq: F, der: 1 - F }; };
  C['chi-no-central'] = v => { need(fin(v.x, v.df, v.lam), 'Rellena todos los campos.'); pos(v.df, 'Los gl deben ser mayores que 0.'); need(v.lam >= 0 && v.x >= 0, 'χ² y λ no pueden ser negativos.');
    const F = U.ncChiCdf(v.x, v.df, v.lam); return { izq: F, der: 1 - F, media: v.df + v.lam }; };

  // Probabilidad y funciones especiales
  const prob = (x, m) => need(x >= 0 && x <= 1, m || 'Las probabilidades deben estar entre 0 y 1.');
  C['probabilidad-sucesos'] = v => { prob(v.pa); prob(v.pb);
    const inter = v.pa * v.pb; return { inter, union: v.pa + v.pb - inter, noa: 1 - v.pa, ninguno: (1 - v.pa) * (1 - v.pb), solouno: v.pa + v.pb - 2 * inter }; };
  C['probabilidad-condicional'] = v => { prob(v.pab); prob(v.pb); need(v.pb > 0, 'P(B) debe ser mayor que 0.'); need(v.pab <= v.pb, 'P(A y B) no puede ser mayor que P(B).');
    return { pcond: v.pab / v.pb }; };
  C['teorema-bayes-diagnostico'] = v => { prob(v.sens); prob(v.esp); prob(v.prev);
    const vp = v.sens * v.prev, fp = (1 - v.esp) * (1 - v.prev), fn = (1 - v.sens) * v.prev, vn = v.esp * (1 - v.prev);
    return { vpp: vp / (vp + fp), vpn: vn / (vn + fn), rvp: v.esp < 1 ? v.sens / (1 - v.esp) : Infinity, rvn: (1 - v.sens) / v.esp, ppos: vp + fp }; };
  C['combinatoria'] = v => { intk(v.n, 'n debe ser un entero no negativo.'); intk(v.k, 'k debe ser un entero no negativo.'); need(v.k <= v.n, 'k no puede superar n.'); need(v.n <= 170, 'n debe ser como máximo 170.');
    return { fact: J.factorial(v.n), comb: Math.round(Math.exp(U.lnchoose(v.n, v.k))), var: Math.round(Math.exp(U.lnfact(v.n) - U.lnfact(v.n - v.k))), varrep: v.n ** v.k, combrep: Math.round(Math.exp(U.lnchoose(v.n + v.k - 1, v.k))) }; };
  C['funcion-gamma'] = v => { need(fin(v.x) && v.x > 0, 'x debe ser mayor que 0.'); return { g: J.gammafn(v.x), lg: J.gammaln(v.x) }; };
  C['gamma-incompleta'] = v => { pos(v.a, 'a debe ser mayor que 0.'); need(v.x >= 0, 'x no puede ser negativo.');
    const P = J.lowRegGamma(v.a, v.x), G = J.gammafn(v.a); return { P, Q: 1 - P, inf: P * G, sup: (1 - P) * G }; };
  C['funcion-beta'] = v => { pos(v.a, 'a debe ser mayor que 0.'); pos(v.b, 'b debe ser mayor que 0.'); return { B: J.betafn(v.a, v.b), lB: J.betaln(v.a, v.b) }; };
  C['beta-incompleta'] = v => { pos(v.a, 'a debe ser mayor que 0.'); pos(v.b, 'b debe ser mayor que 0.'); need(v.x >= 0 && v.x <= 1, 'x debe estar entre 0 y 1.');
    const I = J.ibeta(v.x, v.a, v.b); return { I, inf: I * J.betafn(v.a, v.b), comp: 1 - I }; };
  C['funcion-error'] = v => { need(fin(v.x), 'Escribe un valor.'); return { erf: J.erf(v.x), erfc: J.erfc(v.x) }; };
  C['funcion-error-inversa'] = v => { need(v.y > -1 && v.y < 1, 'El valor debe estar entre −1 y 1.'); const inv = J.erfcinv(1 - v.y); return { inv, z: inv * Math.SQRT2 }; };
})(typeof window !== 'undefined' ? window : globalThis);

(function (root) {
  const C = root.CALC, U = root.CALC_U, J = U.J, { need, fin, mean, vari, sum } = U;
  const lista = (a, min, m) => need(Array.isArray(a) && a.length >= min && a.every(Number.isFinite), m || `Introduce al menos ${min} valores numéricos.`);

  // Descriptivos
  C['descriptivos'] = (v, H) => { lista(v.datos, 3); const a = v.datos, n = a.length, m = mean(a), s = Math.sqrt(vari(a));
    const m2 = sum(a.map(x => (x - m) ** 2)) / n, m3 = sum(a.map(x => (x - m) ** 3)) / n, m4 = sum(a.map(x => (x - m) ** 4)) / n;
    const g1 = m3 / m2 ** 1.5, G1 = n > 2 ? g1 * Math.sqrt(n * (n - 1)) / (n - 2) : NaN;
    const g2 = m4 / m2 ** 2 - 3, G2 = n > 3 ? ((n - 1) / ((n - 2) * (n - 3))) * ((n + 1) * g2 + 6) : NaN;
    return { n, media: m, mediana: U.quantile(a, .5), de: s, var: s * s, ee: s / Math.sqrt(n), min: Math.min(...a), max: Math.max(...a), q1: U.quantile(a, .25), q3: U.quantile(a, .75), asim: G1, curt: G2,
      apa: `M = ${H.f(m)}, DT = ${H.f(s)}, Mdn = ${H.f(U.quantile(a, .5))}, n = ${n}` }; };
  C['puntuaciones-tipicas'] = v => { need(fin(v.x, v.mu, v.sigma), 'Rellena todos los campos.'); need(v.sigma > 0, 'σ debe ser mayor que 0.');
    const z = (v.x - v.mu) / v.sigma; return { z, t: 50 + 10 * z, ci: 100 + 15 * z, ene: Math.min(9, Math.max(1, Math.round(2 * z + 5))), deca: Math.min(10, Math.max(1, Math.round(2 * z + 5.5))), pc: 100 * U.Phi(z) }; };
  C['puntuacion-z'] = v => { need(fin(v.x, v.mu, v.sigma), 'Rellena todos los campos.'); need(v.sigma > 0, 'σ debe ser mayor que 0.');
    const z = (v.x - v.mu) / v.sigma; return { z, pc: 100 * U.Phi(z), sup: 100 * (1 - U.Phi(z)) }; };
  C['coeficiente-variacion'] = v => { need(fin(v.m, v.s), 'Rellena todos los campos.'); need(v.m !== 0, 'La media no puede ser 0.'); need(v.s >= 0, 'La DT no puede ser negativa.');
    return { cv: 100 * v.s / Math.abs(v.m), cvc: Number.isFinite(v.n) && v.n > 1 ? 100 * (1 + 1 / (4 * v.n)) * v.s / Math.abs(v.m) : NaN }; };
  C['valores-atipicos'] = (v, H) => { lista(v.datos, 4); const a = v.datos, m = mean(a), s = Math.sqrt(vari(a)), q1 = U.quantile(a, .25), q3 = U.quantile(a, .75), iqr = q3 - q1;
    const lo = q1 - 1.5 * iqr, hi = q3 + 1.5 * iqr, med = U.quantile(a, .5), mad = 1.4826 * U.quantile(a.map(x => Math.abs(x - med)), .5);
    const tuk = a.filter(x => x < lo || x > hi), zz = a.filter(x => Math.abs((x - m) / s) > 3), mm = mad > 0 ? a.filter(x => Math.abs(x - med) / mad > 3) : [];
    const lst = x => x.length ? x.map(y => H.f(y)).join('; ') : 'Ninguno';
    return { vallas: `[${H.f(lo)}; ${H.f(hi)}]`, tukey: lst(tuk), zeta: lst(zz), mad: lst(mm) }; };

  // Contrastes con medias
  C['p-valor'] = (v, H) => { const d = v.dist; need(fin(v.x), 'Escribe el valor del estadístico.');
    if (d !== 'z') need(v.df1 > 0, 'Los grados de libertad deben ser mayores que 0.'); if (d === 'f') need(v.df2 > 0, 'Los gl del denominador deben ser mayores que 0.');
    if (d === 'chi' || d === 'f') need(v.x >= 0, 'χ² y F no pueden ser negativos.');
    const c = +v.colas; let p, lab;
    if (d === 't') { p = c === 2 ? U.pT(v.x, v.df1) : 1 - J.studentt.cdf(v.x, v.df1); lab = `t(${v.df1}) = ${H.f(v.x)}`; }
    else if (d === 'z') { p = c === 2 ? U.pTwo(v.x) : 1 - U.Phi(v.x); lab = `z = ${H.f(v.x)}`; }
    else if (d === 'chi') { p = 1 - J.chisquare.cdf(v.x, v.df1); lab = `χ²(${v.df1}) = ${H.f(v.x)}`; }
    else { p = 1 - J.centralF.cdf(v.x, v.df1, v.df2); lab = `F(${v.df1}, ${v.df2}) = ${H.f(v.x)}`; }
    return { p, dec: p < .05 ? 'Significativo (α = ,05)' : 'No significativo (α = ,05)', apa: `${lab}, ${H.apaP(p)}` }; };
  C['t-una-muestra'] = (v, H) => { need(fin(v.m, v.s, v.n, v.mu0), 'Rellena todos los campos.'); need(v.s > 0, 'La DT debe ser mayor que 0.'); need(v.n >= 2, 'Se necesitan al menos 2 casos.');
    const se = v.s / Math.sqrt(v.n), t = (v.m - v.mu0) / se, df = v.n - 1, p = U.pT(t, df), tc = J.studentt.inv(.975, df), d = (v.m - v.mu0) / v.s;
    return { t, df, p, d, ic: H.ci(v.m - v.mu0 - tc * se, v.m - v.mu0 + tc * se), apa: `t(${df}) = ${H.f(t)}, ${H.apaP(p)}, d = ${H.f(d)}` }; };
  C['t-independientes'] = (v, H) => { need(fin(v.m1, v.s1, v.n1, v.m2, v.s2, v.n2), 'Rellena todos los campos.'); need(v.s1 > 0 && v.s2 > 0, 'Las DT deben ser mayores que 0.'); need(v.n1 >= 2 && v.n2 >= 2, 'Cada grupo necesita al menos 2 casos.');
    const dif = v.m1 - v.m2, N = v.n1 + v.n2, sp = Math.sqrt(((v.n1 - 1) * v.s1 ** 2 + (v.n2 - 1) * v.s2 ** 2) / (N - 2));
    const seS = sp * Math.sqrt(1 / v.n1 + 1 / v.n2), tS = dif / seS, dfS = N - 2;
    const a = v.s1 ** 2 / v.n1, b = v.s2 ** 2 / v.n2, seW = Math.sqrt(a + b), tW = dif / seW, dfW = (a + b) ** 2 / (a * a / (v.n1 - 1) + b * b / (v.n2 - 1));
    const welch = v.metodo === 'welch', t = welch ? tW : tS, df = welch ? dfW : dfS, se = welch ? seW : seS, p = U.pT(t, df), tc = J.studentt.inv(.975, df);
    return { t, df: H.f(df, welch ? 2 : 0), p, d: dif / sp, ic: H.ci(dif - tc * se, dif + tc * se), apa: `t(${H.f(df, welch ? 2 : 0)}) = ${H.f(t)}, ${H.apaP(p)}, d = ${H.f(dif / sp)}` }; };
  C['t-relacionadas'] = (v, H) => { need(fin(v.md, v.sd, v.n), 'Rellena todos los campos.'); need(v.sd > 0, 'La DT de las diferencias debe ser mayor que 0.'); need(v.n >= 2, 'Se necesitan al menos 2 pares.');
    const se = v.sd / Math.sqrt(v.n), t = v.md / se, df = v.n - 1, p = U.pT(t, df), tc = J.studentt.inv(.975, df), dz = v.md / v.sd;
    return { t, df, p, dz, ic: H.ci(v.md - tc * se, v.md + tc * se), apa: `t(${df}) = ${H.f(t)}, ${H.apaP(p)}, dz = ${H.f(dz)}` }; };
  C['t-datos-brutos'] = (v, H) => { lista(v.g1, 2, 'El grupo 1 necesita al menos 2 valores.'); lista(v.g2, 2, 'El grupo 2 necesita al menos 2 valores.');
    const r = C['t-independientes']({ m1: mean(v.g1), s1: Math.sqrt(vari(v.g1)), n1: v.g1.length, m2: mean(v.g2), s2: Math.sqrt(vari(v.g2)), n2: v.g2.length, metodo: v.metodo }, H);
    return { m1: mean(v.g1), m2: mean(v.g2), s1: Math.sqrt(vari(v.g1)), s2: Math.sqrt(vari(v.g2)), ...r }; };
  C['anova-resumen'] = (v, H) => { lista(v.medias, 2, 'Introduce al menos 2 medias.'); lista(v.des, 2, 'Introduce las DT.'); lista(v.ns, 2, 'Introduce los tamaños.');
    const k = v.medias.length; need(v.des.length === k && v.ns.length === k, 'Debe haber el mismo número de medias, DT y tamaños.'); need(v.ns.every(n => n >= 2), 'Cada grupo necesita al menos 2 casos.');
    const N = sum(v.ns), gm = sum(v.medias.map((m, i) => m * v.ns[i])) / N;
    const ssb = sum(v.medias.map((m, i) => v.ns[i] * (m - gm) ** 2)), ssw = sum(v.des.map((s, i) => (v.ns[i] - 1) * s * s));
    const df1 = k - 1, df2 = N - k, msw = ssw / df2, F = (ssb / df1) / msw, p = 1 - J.centralF.cdf(F, df1, df2), eta = ssb / (ssb + ssw), om = Math.max(0, (ssb - df1 * msw) / (ssb + ssw + msw));
    return { F, gl: `${df1}, ${df2}`, p, eta, omega: om, apa: `F(${df1}, ${df2}) = ${H.f(F)}, ${H.apaP(p)}, η² = ${H.sz(H.f(eta))}, ω² = ${H.sz(H.f(om))}` }; };
  C['comparar-varianzas'] = (v, H) => { need(fin(v.s1, v.n1, v.s2, v.n2), 'Rellena todos los campos.'); need(v.s1 > 0 && v.s2 > 0, 'Las DT deben ser mayores que 0.'); need(v.n1 >= 2 && v.n2 >= 2, 'Cada grupo necesita al menos 2 casos.');
    const F = v.s1 ** 2 / v.s2 ** 2, d1 = v.n1 - 1, d2 = v.n2 - 1, c = J.centralF.cdf(F, d1, d2), p = 2 * Math.min(c, 1 - c);
    return { F, gl: `${d1}, ${d2}`, p, apa: `F(${d1}, ${d2}) = ${H.f(F)}, ${H.apaP(p)}` }; };

  // Proporciones y tablas
  C['z-una-proporcion'] = (v, H) => { need(Number.isInteger(v.x) && Number.isInteger(v.n) && v.n > 0 && v.x >= 0 && v.x <= v.n, 'x y n deben ser enteros con 0 ≤ x ≤ n.'); need(v.p0 > 0 && v.p0 < 1, 'p₀ debe estar entre 0 y 1.');
    const ph = v.x / v.n, z = (ph - v.p0) / Math.sqrt(v.p0 * (1 - v.p0) / v.n), p = U.pTwo(z), h = 2 * Math.asin(Math.sqrt(ph)) - 2 * Math.asin(Math.sqrt(v.p0));
    return { ph, z, p, h, apa: `z = ${H.f(z)}, ${H.apaP(p)}` }; };
  C['z-dos-proporciones'] = (v, H) => { need([v.x1, v.n1, v.x2, v.n2].every(Number.isInteger) && v.n1 > 0 && v.n2 > 0 && v.x1 <= v.n1 && v.x2 <= v.n2 && v.x1 >= 0 && v.x2 >= 0, 'Los valores deben ser enteros con 0 ≤ x ≤ n.');
    const p1 = v.x1 / v.n1, p2 = v.x2 / v.n2, pp = (v.x1 + v.x2) / (v.n1 + v.n2), z = (p1 - p2) / Math.sqrt(pp * (1 - pp) * (1 / v.n1 + 1 / v.n2)), p = U.pTwo(z);
    const se = Math.sqrt(p1 * (1 - p1) / v.n1 + p2 * (1 - p2) / v.n2); return { p1, p2, z, p, ic: H.ci(p1 - p2 - 1.959964 * se, p1 - p2 + 1.959964 * se, 3), apa: `z = ${H.f(z)}, ${H.apaP(p)}` }; };
  C['prueba-binomial'] = (v, H) => { need(Number.isInteger(v.x) && Number.isInteger(v.n) && v.n > 0 && v.x >= 0 && v.x <= v.n, 'x y n deben ser enteros con 0 ≤ x ≤ n.'); need(v.p0 > 0 && v.p0 < 1, 'p₀ debe estar entre 0 y 1.');
    const pk = k => J.binomial.pdf(k, v.n, v.p0), po = pk(v.x); let p2 = 0; for (let k = 0; k <= v.n; k++) { const q = pk(k); if (q <= po * (1 + 1e-7)) p2 += q; }
    const up = 1 - (v.x > 0 ? J.binomial.cdf(v.x - 1, v.n, v.p0) : 0), lo = J.binomial.cdf(v.x, v.n, v.p0);
    return { ph: v.x / v.n, p: Math.min(1, p2), pmay: up, pmen: lo, apa: `prueba binomial exacta, ${H.apaP(Math.min(1, p2))}` }; };
  C['chi-bondad-ajuste'] = (v, H) => { lista(v.obs, 2, 'Introduce al menos 2 frecuencias observadas.'); const o = v.obs, k = o.length, N = sum(o);
    let e; if (Array.isArray(v.esp) && v.esp.length) { need(v.esp.length === k, 'Debe haber tantas proporciones esperadas como categorías.'); const s = sum(v.esp); e = v.esp.map(p => N * p / s); } else e = o.map(() => N / k);
    const chi = sum(o.map((x, i) => (x - e[i]) ** 2 / e[i])), df = k - 1, p = 1 - J.chisquare.cdf(chi, df), w = Math.sqrt(chi / N);
    return { chi, df, p, w, emin: Math.min(...e), apa: `χ²(${df}, N = ${N}) = ${H.f(chi)}, ${H.apaP(p)}, w = ${H.f(w)}` }; };
  const tabla = (m, msg) => { need(Array.isArray(m) && m.length >= 2 && m.every(r => Array.isArray(r) && r.length === m[0].length && r.length >= 2 && r.every(x => Number.isFinite(x) && x >= 0)), msg || 'Introduce una tabla de al menos 2 × 2 con frecuencias no negativas, una fila por línea.'); };
  const chiTabla = m => { const R = m.map(sum), Cc = m[0].map((_, j) => sum(m.map(r => r[j]))), N = sum(R); let chi = 0, emin = Infinity;
    m.forEach((r, i) => r.forEach((x, j) => { const e = R[i] * Cc[j] / N; emin = Math.min(emin, e); chi += (x - e) ** 2 / e; })); return { chi, df: (m.length - 1) * (m[0].length - 1), N, emin, R, Cc }; };
  C['chi-tabla-2x2'] = (v, H) => { const m = [[v.a, v.b], [v.c, v.d]]; tabla(m, 'Las cuatro frecuencias deben ser números no negativos.');
    const t = chiTabla(m), N = t.N, num = Math.abs(v.a * v.d - v.b * v.c), den = t.R[0] * t.R[1] * t.Cc[0] * t.Cc[1];
    const yates = N * Math.max(0, num - N / 2) ** 2 / den, chi = v.yates === 'si' ? yates : t.chi, p = 1 - J.chisquare.cdf(chi, 1), phi = (v.a * v.d - v.b * v.c) / Math.sqrt(den);
    return { chi, p, phi, or: (v.a * v.d) / (v.b * v.c), emin: t.emin, apa: `χ²(1, N = ${N}) = ${H.f(chi)}, ${H.apaP(p)}, φ = ${H.sz(H.f(phi))}` }; };
  C['chi-tabla-rxc'] = (v, H) => { tabla(v.tabla); const t = chiTabla(v.tabla), p = 1 - J.chisquare.cdf(t.chi, t.df), q = Math.min(v.tabla.length, v.tabla[0].length) - 1, V = Math.sqrt(t.chi / (t.N * q));
    return { chi: t.chi, df: t.df, p, V, emin: t.emin, apa: `χ²(${t.df}, N = ${t.N}) = ${H.f(t.chi)}, ${H.apaP(p)}, V = ${H.sz(H.f(V))}` }; };
  C['fisher-exacta-2x2'] = (v, H) => { const m = [v.a, v.b, v.c, v.d]; need(m.every(x => Number.isInteger(x) && x >= 0), 'Las frecuencias deben ser enteros no negativos.');
    const r1 = v.a + v.b, r2 = v.c + v.d, c1 = v.a + v.c, N = r1 + r2, lo = Math.max(0, c1 - r2), hi = Math.min(r1, c1);
    const pr = x => Math.exp(U.lnchoose(r1, x) + U.lnchoose(r2, c1 - x) - U.lnchoose(N, c1)), po = pr(v.a); let p2 = 0, pl = 0, pu = 0;
    for (let x = lo; x <= hi; x++) { const q = pr(x); if (q <= po * (1 + 1e-7)) p2 += q; if (x <= v.a) pl += q; if (x >= v.a) pu += q; }
    return { p: Math.min(1, p2), pl: Math.min(1, pl), pu: Math.min(1, pu), or: (v.a * v.d) / (v.b * v.c), apa: `prueba exacta de Fisher, ${H.apaP(Math.min(1, p2))}` }; };
  C['mcnemar'] = (v, H) => { need(Number.isInteger(v.b) && Number.isInteger(v.c) && v.b >= 0 && v.c >= 0 && v.b + v.c > 0, 'b y c deben ser enteros no negativos y no ambos 0.');
    const n = v.b + v.c, chi = (v.b - v.c) ** 2 / n, chic = (Math.abs(v.b - v.c) - 1) ** 2 / n, p = 1 - J.chisquare.cdf(chi, 1);
    const k = Math.min(v.b, v.c), pex = Math.min(1, 2 * J.binomial.cdf(k, n, .5));
    return { chi, p, chic, pex, apa: `χ²(1) de McNemar = ${H.f(chi)}, ${H.apaP(p)}` }; };

  // No paramétricas
  C['mann-whitney'] = (v, H) => { need(fin(v.u, v.n1, v.n2) && v.n1 > 0 && v.n2 > 0 && v.u >= 0 && v.u <= v.n1 * v.n2, 'Revisa U, n₁ y n₂ (0 ≤ U ≤ n₁·n₂).');
    const mu = v.n1 * v.n2 / 2, s = Math.sqrt(v.n1 * v.n2 * (v.n1 + v.n2 + 1) / 12), z = (v.u - mu) / s, p = U.pTwo(z), r = Math.abs(z) / Math.sqrt(v.n1 + v.n2);
    return { z, p, r, ps: v.u / (v.n1 * v.n2), rbis: 1 - 2 * v.u / (v.n1 * v.n2), apa: `U = ${H.f(v.u, 0)}, z = ${H.f(z)}, ${H.apaP(p)}, r = ${H.sz(H.f(r))}` }; };
  C['wilcoxon'] = (v, H) => { need(fin(v.t, v.n) && v.n > 0 && v.t >= 0 && v.t <= v.n * (v.n + 1) / 2, 'Revisa T y n (0 ≤ T ≤ n(n+1)/2).');
    const mu = v.n * (v.n + 1) / 4, s = Math.sqrt(v.n * (v.n + 1) * (2 * v.n + 1) / 24), z = (v.t - mu) / s, p = U.pTwo(z), r = Math.abs(z) / Math.sqrt(v.n);
    return { z, p, r, rbis: 1 - 4 * v.t / (v.n * (v.n + 1)), apa: `T = ${H.f(v.t, 0)}, z = ${H.f(z)}, ${H.apaP(p)}, r = ${H.sz(H.f(r))}` }; };
  C['kruskal-wallis'] = (v, H) => { need(fin(v.h, v.k, v.n) && v.h >= 0 && v.k >= 2 && v.n > v.k, 'Revisa H, k y N.');
    const df = v.k - 1, p = 1 - J.chisquare.cdf(v.h, df), eps = v.h * (v.n + 1) / (v.n * v.n - 1), eta = (v.h - v.k + 1) / (v.n - v.k);
    return { df, p, eps, eta: Math.max(0, eta), apa: `H(${df}) = ${H.f(v.h)}, ${H.apaP(p)}, ε² = ${H.sz(H.f(eps))}` }; };
})(typeof window !== 'undefined' ? window : globalThis);

(function (root) {
  const C = root.CALC, U = root.CALC_U, J = U.J, { need, fin, mean, vari, sum } = U;
  const rr = (r, m) => need(Number.isFinite(r) && r > -1 && r < 1, m || 'Las correlaciones deben estar entre −1 y 1.');
  const lista = (a, min, m) => need(Array.isArray(a) && a.length >= min && a.every(Number.isFinite), m || `Introduce al menos ${min} valores numéricos.`);
  const r2ok = x => need(Number.isFinite(x) && x >= 0 && x < 1, 'R² debe estar entre 0 y 1.');

  C['correlacion-p-valor'] = (v, H) => { rr(v.r); need(v.n >= 3, 'n debe ser al menos 3.');
    const df = v.n - 2, t = v.r * Math.sqrt(df / (1 - v.r ** 2)), p = U.pT(t, df);
    return { t, df, p, r2: v.r ** 2, apa: `r(${df}) = ${H.sz(H.f(v.r))}, ${H.apaP(p)}` }; };
  C['ic-correlacion'] = (v, H) => { rr(v.r); need(v.n >= 4, 'n debe ser al menos 4.'); const { z, pc } = U.conf(v);
    const zr = Math.atanh(v.r), se = 1 / Math.sqrt(v.n - 3), lo = Math.tanh(zr - z * se), hi = Math.tanh(zr + z * se);
    return { zr, se, lo, hi, apa: `r = ${H.sz(H.f(v.r))}, IC ${pc} % [${H.sz(H.f(lo))}; ${H.sz(H.f(hi))}]` }; };
  C['dos-correlaciones'] = (v, H) => { rr(v.r1); rr(v.r2); need(v.n1 >= 4 && v.n2 >= 4, 'Cada muestra necesita al menos 4 casos.');
    const z1 = Math.atanh(v.r1), z2 = Math.atanh(v.r2), z = (z1 - z2) / Math.sqrt(1 / (v.n1 - 3) + 1 / (v.n2 - 3)), p = U.pTwo(z);
    return { z1, z2, z, p, q: z1 - z2, apa: `z = ${H.f(z)}, ${H.apaP(p)}, q = ${H.f(z1 - z2)}` }; };
  C['correlaciones-dependientes'] = (v, H) => { rr(v.r12); rr(v.r13); rr(v.r23); need(v.n >= 4, 'n debe ser al menos 4.');
    const rm = (v.r12 + v.r13) / 2, psi = v.r23 * (1 - 2 * rm ** 2) - 0.5 * rm ** 2 * (1 - 2 * rm ** 2 - v.r23 ** 2), s = psi / (1 - rm ** 2) ** 2;
    const z = (Math.atanh(v.r12) - Math.atanh(v.r13)) * Math.sqrt(v.n - 3) / Math.sqrt(2 - 2 * s), p = U.pTwo(z);
    return { z, p, dif: v.r12 - v.r13, apa: `Z de Steiger = ${H.f(z)}, ${H.apaP(p)}` }; };
  C['correlacion-parcial'] = (v, H) => { rr(v.rxy); rr(v.rxz); rr(v.ryz); need(v.n >= 5, 'n debe ser al menos 5.');
    const num = v.rxy - v.rxz * v.ryz, pr = num / Math.sqrt((1 - v.rxz ** 2) * (1 - v.ryz ** 2)), sr = num / Math.sqrt(1 - v.ryz ** 2);
    need(Math.abs(pr) < 1, 'Las correlaciones introducidas no son compatibles entre sí.');
    const df = v.n - 3, t = pr * Math.sqrt(df / (1 - pr ** 2)), p = U.pT(t, df);
    return { pr, sr, sr2: sr * sr, t, p, apa: `r parcial(${df}) = ${H.sz(H.f(pr))}, ${H.apaP(p)}` }; };
  C['correlacion-datos'] = (v, H) => { lista(v.x, 3); lista(v.y, 3); need(v.x.length === v.y.length, 'X e Y deben tener el mismo número de valores.');
    const n = v.x.length, r = U.pearson(v.x, v.y), rho = U.pearson(U.ranks(v.x), U.ranks(v.y)), df = n - 2;
    const t = r * Math.sqrt(df / (1 - r * r)), p = U.pT(t, df), ts = rho * Math.sqrt(df / (1 - rho * rho)), ps = U.pT(ts, df);
    return { n, r, p, rho, ps, apa: `r(${df}) = ${H.sz(H.f(r))}, ${H.apaP(p)}` }; };
  C['covarianza-correlacion'] = v => { need(fin(v.sx, v.sy) && v.sx > 0 && v.sy > 0, 'Las DT deben ser mayores que 0.');
    if (v.modo === 'cov') { need(fin(v.val), 'Escribe la covarianza.'); const r = v.val / (v.sx * v.sy); need(Math.abs(r) <= 1, 'La covarianza es incompatible con esas DT (|r| > 1).'); return { r, cov: v.val }; }
    rr(v.val, 'r debe estar entre −1 y 1.'); return { r: v.val, cov: v.val * v.sx * v.sy }; };

  C['regresion-simple-datos'] = (v, H) => { lista(v.x, 3); lista(v.y, 3); need(v.x.length === v.y.length, 'X e Y deben tener el mismo número de valores.');
    const n = v.x.length, mx = mean(v.x), my = mean(v.y); let sxy = 0, sxx = 0; for (let i = 0; i < n; i++) { sxy += (v.x[i] - mx) * (v.y[i] - my); sxx += (v.x[i] - mx) ** 2; }
    need(sxx > 0, 'X debe tener variabilidad.'); const b1 = sxy / sxx, b0 = my - b1 * mx, res = v.y.map((y, i) => y - b0 - b1 * v.x[i]), sse = sum(res.map(e => e * e)), sst = sum(v.y.map(y => (y - my) ** 2));
    const se = Math.sqrt(sse / (n - 2)), seb = se / Math.sqrt(sxx), t = b1 / seb, p = U.pT(t, n - 2), R2 = 1 - sse / sst;
    return { b0, b1, R2, see: se, t, p, apa: `b = ${H.f(b1)}, t(${n - 2}) = ${H.f(t)}, ${H.apaP(p)}, R² = ${H.sz(H.f(R2))}` }; };
  C['r2-ajustado'] = v => { r2ok(v.r2); need(v.n > v.k + 1 && v.k >= 1, 'n debe ser mayor que k + 1.');
    const adj = 1 - (1 - v.r2) * (v.n - 1) / (v.n - v.k - 1); const op = 1 - ((v.n - 3) / (v.n - v.k - 1)) * (1 - v.r2) * (1 + 2 * (1 - v.r2) / (v.n - v.k + 1));
    return { adj, op: Math.max(0, op), contr: v.r2 - adj }; };
  C['f-r2'] = (v, H) => { r2ok(v.r2); need(v.n > v.k + 1 && v.k >= 1, 'n debe ser mayor que k + 1.');
    const df2 = v.n - v.k - 1, F = (v.r2 / v.k) / ((1 - v.r2) / df2), p = 1 - J.centralF.cdf(F, v.k, df2), f2 = v.r2 / (1 - v.r2);
    return { F, gl: `${v.k}, ${df2}`, p, f2, apa: `R² = ${H.sz(H.f(v.r2))}, F(${v.k}, ${df2}) = ${H.f(F)}, ${H.apaP(p)}` }; };
  C['cambio-r2'] = (v, H) => { r2ok(v.r2a); r2ok(v.r2b); need(v.r2b >= v.r2a, 'R² del modelo ampliado debe ser mayor o igual que el del modelo reducido.');
    need(v.ka >= 1 && v.kb >= 1 && v.n > v.ka + v.kb + 1, 'Revisa n y el número de predictores.');
    const df2 = v.n - v.ka - v.kb - 1, dR = v.r2b - v.r2a, F = (dR / v.kb) / ((1 - v.r2b) / df2), p = 1 - J.centralF.cdf(F, v.kb, df2), f2 = dR / (1 - v.r2b);
    return { dR, F, gl: `${v.kb}, ${df2}`, p, f2, apa: `ΔR² = ${H.sz(H.f(dR, 3))}, F(${v.kb}, ${df2}) = ${H.f(F)}, ${H.apaP(p)}, f² = ${H.f(f2)}` }; };
  C['ic-r2'] = (v, H) => { r2ok(v.r2); need(v.n > v.k + 3 && v.k >= 1, 'n debe ser mayor que k + 3.'); const { z, pc } = U.conf(v);
    const se = Math.sqrt(4 * v.r2 * (1 - v.r2) ** 2 * (v.n - v.k - 1) ** 2 / ((v.n * v.n - 1) * (v.n + 3)));
    const lo = Math.max(0, v.r2 - z * se), hi = Math.min(1, v.r2 + z * se);
    return { se, lo, hi, apa: `R² = ${H.sz(H.f(v.r2))}, IC ${pc} % [${H.sz(H.f(lo))}; ${H.sz(H.f(hi))}]` }; };
  C['ic-coeficiente'] = (v, H) => { need(fin(v.b, v.se, v.n, v.k) && v.se > 0, 'Rellena todos los campos (EE > 0).'); need(v.n > v.k + 1, 'n debe ser mayor que k + 1.');
    const df = v.n - v.k - 1, tc = J.studentt.inv(1 - (1 - v.conf) / 2, df), t = v.b / v.se, p = U.pT(t, df), pc = Math.round(v.conf * 100);
    return { t, df, p, lo: v.b - tc * v.se, hi: v.b + tc * v.se, apa: `b = ${H.f(v.b)}, EE = ${H.f(v.se)}, t(${df}) = ${H.f(t)}, ${H.apaP(p)}, IC ${pc} % [${H.f(v.b - tc * v.se)}; ${H.f(v.b + tc * v.se)}]` }; };
  C['comparar-pendientes'] = (v, H) => { need(fin(v.b1, v.se1, v.n1, v.b2, v.se2, v.n2) && v.se1 > 0 && v.se2 > 0, 'Rellena todos los campos (EE > 0).'); need(v.n1 > 2 && v.n2 > 2, 'Cada muestra necesita más de 2 casos.');
    const df = v.n1 + v.n2 - 4, t = (v.b1 - v.b2) / Math.sqrt(v.se1 ** 2 + v.se2 ** 2), p = U.pT(t, df);
    return { dif: v.b1 - v.b2, t, df, p, apa: `t(${df}) = ${H.f(t)}, ${H.apaP(p)}` }; };
  C['ic-prediccion'] = (v, H) => { need(fin(v.b0, v.b1, v.see, v.n, v.mx, v.sx, v.x0) && v.see > 0 && v.sx > 0 && v.n > 2, 'Rellena todos los campos (EE y DT > 0, n > 2).');
    const yh = v.b0 + v.b1 * v.x0, df = v.n - 2, tc = J.studentt.inv(1 - (1 - v.conf) / 2, df), h = 1 / v.n + (v.x0 - v.mx) ** 2 / ((v.n - 1) * v.sx ** 2);
    const sm = v.see * Math.sqrt(h), sp = v.see * Math.sqrt(1 + h);
    return { yh, icm: H.ci(yh - tc * sm, yh + tc * sm), icp: H.ci(yh - tc * sp, yh + tc * sp) }; };
  C['beta-estandarizado'] = v => { need(fin(v.b, v.sx, v.sy) && v.sx > 0 && v.sy > 0, 'Rellena todos los campos (DT > 0).');
    return { beta: v.b * v.sx / v.sy, bdesde: v.b }; };

  // Mediación
  C['efecto-indirecto-sobel'] = (v, H) => { need(fin(v.a, v.sa, v.b, v.sb) && v.sa > 0 && v.sb > 0, 'Rellena todos los campos (EE > 0).');
    const ab = v.a * v.b, sS = Math.sqrt(v.b ** 2 * v.sa ** 2 + v.a ** 2 * v.sb ** 2), sA = Math.sqrt(sS ** 2 + v.sa ** 2 * v.sb ** 2), sG = Math.sqrt(Math.max(1e-300, sS ** 2 - v.sa ** 2 * v.sb ** 2));
    return { ab, se: sS, z: ab / sS, p: U.pTwo(ab / sS), zA: ab / sA, zG: ab / sG, apa: `ab = ${H.f(ab, 3)}, z de Sobel = ${H.f(ab / sS)}, ${H.apaP(U.pTwo(ab / sS))}` }; };
  C['efecto-indirecto-montecarlo'] = (v, H) => { need(fin(v.a, v.sa, v.b, v.sb) && v.sa > 0 && v.sb > 0, 'Rellena todos los campos (EE > 0).');
    const R = Math.min(200000, Math.max(1000, Math.round(v.reps || 20000))); let s = 20260101;
    const rnd = () => { s = (s * 1664525 + 1013904223) >>> 0; return (s + 0.5) / 4294967296; };
    const norm = () => Math.sqrt(-2 * Math.log(rnd())) * Math.cos(2 * Math.PI * rnd());
    const ab = new Float64Array(R); for (let i = 0; i < R; i++) ab[i] = (v.a + v.sa * norm()) * (v.b + v.sb * norm()); ab.sort();
    const al = 1 - v.conf, lo = ab[Math.floor(R * al / 2)], hi = ab[Math.ceil(R * (1 - al / 2)) - 1], pc = Math.round(v.conf * 100);
    return { ab: v.a * v.b, lo, hi, sig: (lo > 0 || hi < 0) ? 'El IC no incluye 0' : 'El IC incluye 0', apa: `ab = ${H.f(v.a * v.b, 3)}, IC ${pc} % Monte Carlo [${H.f(lo, 3)}; ${H.f(hi, 3)}]` }; };
  C['efecto-total-mediacion'] = v => { need(fin(v.a, v.b, v.cp), 'Rellena todos los campos.'); const ab = v.a * v.b, c = v.cp + ab;
    return { ab, c, pm: c !== 0 ? ab / c : NaN, ratio: v.cp !== 0 ? ab / v.cp : NaN }; };
})(typeof window !== 'undefined' ? window : globalThis);

(function (root) {
  const C = root.CALC, U = root.CALC_U, J = U.J, { need, fin } = U;
  const prop = x => need(Number.isFinite(x) && x >= 0 && x <= 1, 'Las proporciones deben estar entre 0 y 1.');

  // Tamaño del efecto
  C['d-cohen'] = (v, H) => { need(fin(v.m1, v.s1, v.n1, v.m2, v.s2, v.n2), 'Rellena todos los campos.'); need(v.s1 > 0 && v.s2 > 0, 'Las DT deben ser mayores que 0.'); need(v.n1 >= 2 && v.n2 >= 2, 'Cada grupo necesita al menos 2 casos.');
    const N = v.n1 + v.n2, sp = Math.sqrt(((v.n1 - 1) * v.s1 ** 2 + (v.n2 - 1) * v.s2 ** 2) / (N - 2)), d = (v.m1 - v.m2) / sp, Jc = 1 - 3 / (4 * N - 9), g = Jc * d;
    const se = Jc * Math.sqrt(N / (v.n1 * v.n2) + d * d / (2 * N)), { z, pc } = U.conf(v), t = (v.m1 - v.m2) / (sp * Math.sqrt(1 / v.n1 + 1 / v.n2)), df = N - 2, p = U.pT(t, df);
    return { d, g, ic: H.ci(g - z * se, g + z * se), t: `${H.f(t)} (${df})`, p, apa: `t(${df}) = ${H.f(t)}, ${H.apaP(p)}, g = ${H.f(g)}, IC ${pc} % [${H.f(g - z * se)}; ${H.f(g + z * se)}]` }; };
  C['d-ic-exacto'] = (v, H) => { need(fin(v.d, v.n1, v.n2) && v.n1 >= 2 && v.n2 >= 2, 'Rellena todos los campos (n ≥ 2 por grupo).');
    const df = v.n1 + v.n2 - 2, k = Math.sqrt(1 / v.n1 + 1 / v.n2), t = v.d / k, a = 1 - v.conf, pc = Math.round(v.conf * 100);
    const span = Math.abs(t) + 40, lo = U.bisect(n => U.nctCdf(t, df, n) - (1 - a / 2), t - span, t + span), hi = U.bisect(n => U.nctCdf(t, df, n) - a / 2, t - span, t + span);
    return { lo: -(-lo) * k, hi: hi * k, t, apa: `d = ${H.f(v.d)}, IC ${pc} % [${H.f(lo * k)}; ${H.f(hi * k)}]` }; };
  C['d-desde-t'] = (v, H) => { need(fin(v.t, v.n1, v.n2) && v.n1 >= 2 && v.n2 >= 2, 'Rellena todos los campos (n ≥ 2 por grupo).');
    const d = v.t * Math.sqrt(1 / v.n1 + 1 / v.n2), N = v.n1 + v.n2, g = d * (1 - 3 / (4 * N - 9)), r = v.t / Math.sqrt(v.t ** 2 + N - 2);
    return { d, g, r, apa: `d = ${H.f(d)}` }; };
  C['d-relacionadas'] = (v, H) => { need(fin(v.m1, v.m2, v.s1, v.s2, v.r, v.n) && v.s1 > 0 && v.s2 > 0 && Math.abs(v.r) < 1 && v.n >= 2, 'Rellena todos los campos (DT > 0, |r| < 1).');
    const dif = v.m2 - v.m1, sdif = Math.sqrt(v.s1 ** 2 + v.s2 ** 2 - 2 * v.r * v.s1 * v.s2), dz = dif / sdif, dav = dif / ((v.s1 + v.s2) / 2), drm = dz * Math.sqrt(2 * (1 - v.r)), t = dz * Math.sqrt(v.n);
    return { dz, dav, drm, sdif, t, apa: `dz = ${H.f(dz)}, dav = ${H.f(dav)}` }; };
  C['delta-glass'] = (v, H) => { need(fin(v.me, v.mc, v.sc, v.nc) && v.sc > 0, 'Rellena todos los campos (DT > 0).'); const d = (v.me - v.mc) / v.sc;
    const se = Number.isFinite(v.ne) && v.ne > 1 && v.nc > 1 ? Math.sqrt((v.ne + v.nc) / (v.ne * v.nc) + d * d / (2 * (v.nc - 1))) : NaN;
    return { delta: d, ic: Number.isFinite(se) ? H.ci(d - 1.959964 * se, d + 1.959964 * se) : '—', apa: `Δ de Glass = ${H.f(d)}` }; };
  C['eta-desde-f'] = (v, H) => { need(fin(v.F, v.df1, v.df2) && v.F >= 0 && v.df1 > 0 && v.df2 > 0, 'Rellena todos los campos.');
    const ep = v.F * v.df1 / (v.F * v.df1 + v.df2), f = Math.sqrt(ep / (1 - ep)); return { ep, f, apa: `η²p = ${H.sz(H.f(ep))}` }; };
  C['omega-desde-f'] = (v, H) => { need(fin(v.F, v.df1, v.df2) && v.F >= 0 && v.df1 > 0 && v.df2 > 0, 'Rellena todos los campos.');
    const N = v.df1 + v.df2 + 1, om = Math.max(0, v.df1 * (v.F - 1) / (v.df1 * v.F + v.df2 + 1)), eps = Math.max(0, v.df1 * (v.F - 1) / (v.df1 * v.F + v.df2));
    return { om, eps, N, apa: `ω² = ${H.sz(H.f(om))}` }; };
  C['f-cohen-eta'] = v => { if (v.modo === 'eta') { prop(v.val); need(v.val < 1, 'η² debe ser menor que 1.'); return { eta: v.val, f: Math.sqrt(v.val / (1 - v.val)), f2: v.val / (1 - v.val) }; }
    need(fin(v.val) && v.val >= 0, 'f debe ser mayor o igual que 0.'); const e = v.val ** 2 / (1 + v.val ** 2); return { eta: e, f: v.val, f2: v.val ** 2 }; };
  C['conversion-efectos'] = v => { need(fin(v.valor), 'Escribe un valor.'); let d; const x = v.valor;
    if (v.tipo === 'd') d = x; else if (v.tipo === 'r') { need(Math.abs(x) < 1, 'r debe estar entre −1 y 1.'); d = 2 * x / Math.sqrt(1 - x * x); }
    else if (v.tipo === 'or') { need(x > 0, 'La odds ratio debe ser mayor que 0.'); d = Math.log(x) * Math.sqrt(3) / Math.PI; }
    else { need(x >= 0 && x < 1, 'η² debe estar entre 0 y 1.'); d = 2 * Math.sqrt(x) / Math.sqrt(1 - x); }
    const r = d / Math.sqrt(d * d + 4); return { d, r, eta: r * r, f: Math.abs(d) / 2, or: Math.exp(d * Math.PI / Math.sqrt(3)) }; };
  C['riesgo-odds-2x2'] = (v, H) => { const m = [v.a, v.b, v.c, v.d]; need(m.every(x => Number.isFinite(x) && x >= 0) && v.a + v.b > 0 && v.c + v.d > 0, 'Las frecuencias deben ser no negativas.');
    const cc = m.some(x => x === 0) ? 0.5 : 0, [a, b, c, d] = m.map(x => x + cc), { z, pc } = U.conf(v);
    const or = a * d / (b * c), seo = Math.sqrt(1 / a + 1 / b + 1 / c + 1 / d), r1 = a / (a + b), r2 = c / (c + d), rr = r1 / r2, ser = Math.sqrt(1 / a - 1 / (a + b) + 1 / c - 1 / (c + d));
    return { r1, r2, rr, icrr: H.ci(rr * Math.exp(-z * ser), rr * Math.exp(z * ser)), or, icor: H.ci(or * Math.exp(-z * seo), or * Math.exp(z * seo)), dr: r1 - r2,
      apa: `OR = ${H.f(or)}, IC ${pc} % [${H.f(or * Math.exp(-z * seo))}; ${H.f(or * Math.exp(z * seo))}]` }; };
  C['nnt'] = (v, H) => { prop(v.pc); prop(v.pt); need(v.pc !== v.pt, 'Las proporciones deben ser distintas.');
    const rar = v.pc - v.pt, nnt = 1 / Math.abs(rar), rrr = v.pc > 0 ? rar / v.pc : NaN;
    return { rar, nnt: Math.ceil(nnt), nntx: nnt, rrr, tipo: rar > 0 ? 'NNT (beneficio)' : 'NND (daño)' }; };
  C['phi-v-cramer'] = (v, H) => { need(fin(v.chi, v.n, v.f, v.c) && v.chi >= 0 && v.n > 0 && v.f >= 2 && v.c >= 2, 'Rellena todos los campos.');
    const q = Math.min(v.f, v.c) - 1, V = Math.sqrt(v.chi / (v.n * q)), w = Math.sqrt(v.chi / v.n); return { V, w, phi: q === 1 ? w : NaN, apa: `V de Cramér = ${H.sz(H.f(V))}` }; };
  C['w-cohen'] = v => { if (v.modo === 'chi') { need(fin(v.chi, v.n) && v.chi >= 0 && v.n > 0, 'Rellena χ² y N.'); return { w: Math.sqrt(v.chi / v.n) }; }
    need(Array.isArray(v.p0) && Array.isArray(v.p1) && v.p0.length === v.p1.length && v.p0.length >= 2, 'Introduce dos listas de proporciones de igual longitud.');
    const s0 = U.sum(v.p0), s1 = U.sum(v.p1); return { w: Math.sqrt(U.sum(v.p0.map((p, i) => (v.p1[i] / s1 - p / s0) ** 2 / (p / s0)))) }; };
  C['h-cohen'] = (v, H) => { prop(v.p1); prop(v.p2); const h = 2 * Math.asin(Math.sqrt(v.p1)) - 2 * Math.asin(Math.sqrt(v.p2)); return { h, ah: Math.abs(h), apa: `h = ${H.f(h)}` }; };
  C['probabilidad-superioridad'] = v => { need(fin(v.d), 'Escribe d.'); const ps = U.Phi(v.d / Math.SQRT2), u3 = U.Phi(v.d), ovl = 2 * U.Phi(-Math.abs(v.d) / 2);
    return { ps, u3, ovl, nntk: ps > 0.5 ? 1 / (2 * ps - 1) : NaN }; };

  // Estimación
  C['ic-media'] = (v, H) => { need(fin(v.m, v.s, v.n) && v.s > 0 && v.n >= 2, 'Rellena todos los campos (DT > 0, n ≥ 2).');
    const se = v.s / Math.sqrt(v.n), tc = J.studentt.inv(1 - (1 - v.conf) / 2, v.n - 1), pc = Math.round(v.conf * 100);
    return { se, tc, lo: v.m - tc * se, hi: v.m + tc * se, apa: `M = ${H.f(v.m)}, DT = ${H.f(v.s)}, IC ${pc} % [${H.f(v.m - tc * se)}; ${H.f(v.m + tc * se)}]` }; };
  C['ic-media-sigma'] = (v, H) => { need(fin(v.m, v.sigma, v.n) && v.sigma > 0 && v.n >= 1, 'Rellena todos los campos (σ > 0).'); const { z, pc } = U.conf(v), se = v.sigma / Math.sqrt(v.n);
    return { se, z, lo: v.m - z * se, hi: v.m + z * se, apa: `M = ${H.f(v.m)}, IC ${pc} % [${H.f(v.m - z * se)}; ${H.f(v.m + z * se)}]` }; };
  C['ic-proporcion'] = (v, H) => { need(Number.isInteger(v.x) && Number.isInteger(v.n) && v.n > 0 && v.x >= 0 && v.x <= v.n, 'x y n deben ser enteros con 0 ≤ x ≤ n.');
    const { z, c, pc } = U.conf(v), p = v.x / v.n, a = 1 - c, den = 1 + z * z / v.n, cen = (p + z * z / (2 * v.n)) / den, hw = z * Math.sqrt(p * (1 - p) / v.n + z * z / (4 * v.n * v.n)) / den;
    const cpl = v.x === 0 ? 0 : J.beta.inv(a / 2, v.x, v.n - v.x + 1), cpu = v.x === v.n ? 1 : J.beta.inv(1 - a / 2, v.x + 1, v.n - v.x), sw = Math.sqrt(p * (1 - p) / v.n);
    return { p, wilson: H.ci(cen - hw, cen + hw, 3), cp: H.ci(cpl, cpu, 3), wald: H.ci(Math.max(0, p - z * sw), Math.min(1, p + z * sw), 3), apa: `${H.f(100 * p, 1)} %, IC ${pc} % [${H.f(100 * (cen - hw), 1)}; ${H.f(100 * (cen + hw), 1)}]` }; };
  C['ic-diferencia-medias'] = (v, H) => { need(fin(v.m1, v.s1, v.n1, v.m2, v.s2, v.n2) && v.s1 > 0 && v.s2 > 0 && v.n1 >= 2 && v.n2 >= 2, 'Rellena todos los campos.');
    const a = v.s1 ** 2 / v.n1, b = v.s2 ** 2 / v.n2, se = Math.sqrt(a + b), df = (a + b) ** 2 / (a * a / (v.n1 - 1) + b * b / (v.n2 - 1)), tc = J.studentt.inv(1 - (1 - v.conf) / 2, df), dif = v.m1 - v.m2, pc = Math.round(v.conf * 100);
    return { dif, se, df, lo: dif - tc * se, hi: dif + tc * se, apa: `diferencia = ${H.f(dif)}, IC ${pc} % [${H.f(dif - tc * se)}; ${H.f(dif + tc * se)}]` }; };
  C['ic-diferencia-proporciones'] = (v, H) => { need([v.x1, v.n1, v.x2, v.n2].every(Number.isInteger) && v.n1 > 0 && v.n2 > 0 && v.x1 >= 0 && v.x2 >= 0 && v.x1 <= v.n1 && v.x2 <= v.n2, 'Los valores deben ser enteros con 0 ≤ x ≤ n.');
    const { z, pc } = U.conf(v), p1 = v.x1 / v.n1, p2 = v.x2 / v.n2, dif = p1 - p2;
    const wil = (x, n) => { const p = x / n, den = 1 + z * z / n, c = (p + z * z / (2 * n)) / den, h = z * Math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / den; return [c - h, c + h]; };
    const [l1, u1] = wil(v.x1, v.n1), [l2, u2] = wil(v.x2, v.n2), lo = dif - Math.sqrt((p1 - l1) ** 2 + (u2 - p2) ** 2), hi = dif + Math.sqrt((u1 - p1) ** 2 + (p2 - l2) ** 2);
    const se = Math.sqrt(p1 * (1 - p1) / v.n1 + p2 * (1 - p2) / v.n2);
    return { dif, newc: H.ci(lo, hi, 3), wald: H.ci(dif - z * se, dif + z * se, 3), apa: `diferencia = ${H.f(100 * dif, 1)} puntos, IC ${pc} % [${H.f(100 * lo, 1)}; ${H.f(100 * hi, 1)}]` }; };
  C['ic-varianza'] = (v, H) => { need(fin(v.s, v.n) && v.s > 0 && v.n >= 2, 'Rellena todos los campos (DT > 0, n ≥ 2).');
    const df = v.n - 1, a = 1 - v.conf, s2 = v.s ** 2, lo = df * s2 / J.chisquare.inv(1 - a / 2, df), hi = df * s2 / J.chisquare.inv(a / 2, df);
    return { s2, var: H.ci(lo, hi), de: H.ci(Math.sqrt(lo), Math.sqrt(hi)) }; };
  C['ic-poisson'] = (v, H) => { need(Number.isInteger(v.x) && v.x >= 0, 'El recuento debe ser un entero no negativo.'); const a = 1 - v.conf, t = Number.isFinite(v.t) && v.t > 0 ? v.t : 1;
    const lo = v.x === 0 ? 0 : J.chisquare.inv(a / 2, 2 * v.x) / 2, hi = J.chisquare.inv(1 - a / 2, 2 * v.x + 2) / 2;
    return { tasa: v.x / t, lo: lo / t, hi: hi / t }; };
})(typeof window !== 'undefined' ? window : globalThis);

(function (root) {
  const C = root.CALC, U = root.CALC_U, J = U.J, { need, fin, sum, mean } = U;
  const lista = (a, min, m) => need(Array.isArray(a) && a.length >= min && a.every(Number.isFinite), m || `Introduce al menos ${min} valores numéricos.`);
  const rel = (x, m) => need(Number.isFinite(x) && x > 0 && x < 1, m || 'La fiabilidad debe estar entre 0 y 1.');

  // Psicometría
  C['alfa-cronbach'] = (v, H) => { lista(v.vars, 2, 'Introduce la varianza de cada ítem (al menos 2).'); need(v.vt > 0, 'La varianza total debe ser mayor que 0.');
    const k = v.vars.length, a = (k / (k - 1)) * (1 - sum(v.vars) / v.vt); return { k, alfa: a, apa: `α = ${H.sz(H.f(a))}` }; };
  C['alfa-estandarizado'] = (v, H) => { need(Number.isInteger(v.k) && v.k >= 2, 'k debe ser un entero mayor o igual que 2.'); need(v.rm > -1 && v.rm < 1, 'La correlación media debe estar entre −1 y 1.');
    const a = v.k * v.rm / (1 + (v.k - 1) * v.rm); return { alfa: a, apa: `α estandarizado = ${H.sz(H.f(a))}` }; };
  C['ic-alfa-feldt'] = (v, H) => { need(fin(v.alfa, v.n, v.k) && v.alfa < 1 && v.n > 2 && v.k >= 2, 'Rellena todos los campos (α < 1, n > 2, k ≥ 2).');
    const d1 = v.n - 1, d2 = (v.n - 1) * (v.k - 1), a = 1 - v.conf, lo = 1 - (1 - v.alfa) * J.centralF.inv(1 - a / 2, d1, d2), hi = 1 - (1 - v.alfa) * J.centralF.inv(a / 2, d1, d2), pc = Math.round(v.conf * 100);
    return { lo, hi, apa: `α = ${H.sz(H.f(v.alfa))}, IC ${pc} % [${H.sz(H.f(lo))}; ${H.sz(H.f(hi))}]` }; };
  C['omega-mcdonald'] = (v, H) => { lista(v.cargas, 2, 'Introduce al menos 2 cargas factoriales estandarizadas.'); need(v.cargas.every(l => Math.abs(l) < 1), 'Las cargas estandarizadas deben estar entre −1 y 1.');
    const s = sum(v.cargas), e = sum(v.cargas.map(l => 1 - l * l)), om = s * s / (s * s + e), H2 = 1 / (1 + 1 / sum(v.cargas.map(l => l * l / (1 - l * l))));
    return { om, hcoef: H2, ave: mean(v.cargas.map(l => l * l)), apa: `ω = ${H.sz(H.f(om))}` }; };
  C['spearman-brown'] = v => { rel(v.r); need(v.factor > 0, 'El factor de alargamiento debe ser mayor que 0.');
    const rn = v.factor * v.r / (1 + (v.factor - 1) * v.r); let kk = NaN; if (Number.isFinite(v.obj) && v.obj > 0 && v.obj < 1) kk = v.obj * (1 - v.r) / (v.r * (1 - v.obj));
    return { rn, kk, items: Number.isFinite(v.items) ? Math.ceil(v.items * kk) : NaN }; };
  C['error-tipico-medida'] = (v, H) => { need(fin(v.s, v.x) && v.s > 0, 'Rellena la DT y la puntuación.'); rel(v.rxx); const { z, pc } = U.conf(v);
    const sem = v.s * Math.sqrt(1 - v.rxx), see = v.s * Math.sqrt(v.rxx * (1 - v.rxx)), mu = Number.isFinite(v.mu) ? v.mu : v.x, tv = v.rxx * v.x + (1 - v.rxx) * mu;
    return { sem, ico: H.ci(v.x - z * sem, v.x + z * sem), tv, icv: H.ci(tv - z * see, tv + z * see), sed: sem * Math.SQRT2 }; };
  C['cambio-fiable'] = v => { need(fin(v.x1, v.x2, v.s) && v.s > 0, 'Rellena todos los campos (DT > 0).'); rel(v.rxx);
    const sdif = Math.SQRT2 * v.s * Math.sqrt(1 - v.rxx), rci = (v.x2 - v.x1) / sdif;
    return { sdif, rci, minimo: 1.96 * sdif, dec: Math.abs(rci) > 1.96 ? (rci > 0 ? 'Aumento fiable' : 'Descenso fiable') : 'Cambio no fiable' }; };
  C['correccion-atenuacion'] = v => { need(fin(v.r) && Math.abs(v.r) <= 1, 'r debe estar entre −1 y 1.'); rel(v.rxx, 'Las fiabilidades deben estar entre 0 y 1.'); rel(v.ryy, 'Las fiabilidades deben estar entre 0 y 1.');
    const rc = v.r / Math.sqrt(v.rxx * v.ryy); return { rc, rx: v.r / Math.sqrt(v.rxx), max: Math.sqrt(v.rxx * v.ryy) }; };
  C['kappa-cohen'] = (v, H) => { const m = v.tabla; need(Array.isArray(m) && m.length >= 2 && m.every(r => Array.isArray(r) && r.length === m.length && r.every(x => Number.isFinite(x) && x >= 0)), 'Introduce una tabla cuadrada de acuerdo (una fila por línea).');
    const k = m.length, N = sum(m.map(sum)); need(N > 0, 'La tabla está vacía.'); let po = 0, pe = 0;
    for (let i = 0; i < k; i++) { po += m[i][i]; pe += sum(m[i]) * sum(m.map(r => r[i])); } po /= N; pe /= N * N;
    const kap = (po - pe) / (1 - pe), se = Math.sqrt(po * (1 - po) / (N * (1 - pe) ** 2));
    let wo = 0, we = 0; for (let i = 0; i < k; i++) for (let j = 0; j < k; j++) { const w = ((i - j) / (k - 1)) ** 2; wo += w * m[i][j] / N; we += w * sum(m[i]) * sum(m.map(r => r[j])) / (N * N); }
    return { po, pe, kap, ic: H.ci(kap - 1.959964 * se, kap + 1.959964 * se), kw: k > 2 ? 1 - wo / we : NaN, apa: `κ = ${H.sz(H.f(kap))}, IC 95 % [${H.sz(H.f(kap - 1.959964 * se))}; ${H.sz(H.f(kap + 1.959964 * se))}]` }; };
  C['cci-anova'] = (v, H) => { need(fin(v.msb, v.msw, v.k) && v.msb > 0 && v.msw >= 0 && v.k >= 2, 'Rellena todos los campos (k ≥ 2).');
    const i1 = (v.msb - v.msw) / (v.msb + (v.k - 1) * v.msw), ik = (v.msb - v.msw) / v.msb; return { i1, ik, apa: `CCI(1,1) = ${H.sz(H.f(i1))}` }; };

  // Potencia y tamaño muestral
  const pw = x => need(x > 0 && x < 1, 'La potencia debe estar entre 0 y 1.');
  const al = x => need(x > 0 && x < 1, 'α debe estar entre 0 y 1.');
  const dd = x => need(Number.isFinite(x) && x !== 0, 'El tamaño del efecto debe ser distinto de 0.');
  C['tamano-muestral-t'] = (v, H) => { dd(v.d); al(+v.alfa); pw(+v.potencia); const a = +v.alfa, c = +v.colas, target = +v.potencia;
    const n = U.nMin(n => U.powT(Math.abs(v.d), n, a, c, 'ind'), target);
    return { n, N: 2 * n, real: U.powT(Math.abs(v.d), n, a, c, 'ind'), apa: `Para detectar d = ${H.f(Math.abs(v.d))} con α = ${H.sz(H.f(a))} (${c === 2 ? 'bilateral' : 'unilateral'}) y potencia de ${H.sz(H.f(target))} se necesitan ${n} participantes por grupo (N = ${2 * n}).` }; };
  C['potencia-t-independientes'] = v => { dd(v.d); need(v.n >= 2, 'n por grupo debe ser al menos 2.'); const p = U.powT(Math.abs(v.d), v.n, +v.alfa, +v.colas, 'ind'); return { pot: p, beta: 1 - p }; };
  C['tamano-muestral-t-relacionadas'] = (v, H) => { dd(v.d); const n = U.nMin(n => U.powT(Math.abs(v.d), n, +v.alfa, +v.colas, 'una'), +v.potencia);
    return { n, real: U.powT(Math.abs(v.d), n, +v.alfa, +v.colas, 'una'), apa: `Se necesitan ${n} participantes (pares) para detectar dz = ${H.f(Math.abs(v.d))} con potencia de ${H.sz(H.f(+v.potencia))}.` }; };
  C['potencia-t-relacionadas'] = v => { dd(v.d); need(v.n >= 2, 'n debe ser al menos 2.'); const p = U.powT(Math.abs(v.d), v.n, +v.alfa, +v.colas, 'una'); return { pot: p, beta: 1 - p }; };
  const powR = (r, n, a, c) => { const zr = Math.abs(Math.atanh(r)) * Math.sqrt(n - 3), zc = U.zq(1 - a / c); return c === 2 ? U.Phi(zr - zc) + U.Phi(-zr - zc) : U.Phi(zr - zc); };
  C['tamano-muestral-correlacion'] = (v, H) => { need(Math.abs(v.r) > 0 && Math.abs(v.r) < 1, 'r debe estar entre −1 y 1 y ser distinta de 0.');
    const n = U.nMin(n => powR(v.r, n, +v.alfa, +v.colas), +v.potencia, 4); return { n, real: powR(v.r, n, +v.alfa, +v.colas), apa: `Se necesitan ${n} participantes para detectar r = ${H.sz(H.f(v.r))} con potencia de ${H.sz(H.f(+v.potencia))}.` }; };
  C['potencia-correlacion'] = v => { need(Math.abs(v.r) > 0 && Math.abs(v.r) < 1, 'r debe estar entre −1 y 1.'); need(v.n >= 4, 'n debe ser al menos 4.'); const p = powR(v.r, v.n, +v.alfa, +v.colas); return { pot: p, beta: 1 - p }; };
  const powA = (f, k, Ntot, a) => powF(f * f * Ntot, k - 1, Ntot - k, a); const powF = (l, d1, d2, a) => U.powF(l, d1, d2, a);
  C['tamano-muestral-anova'] = (v, H) => { need(v.f > 0, 'f debe ser mayor que 0.'); need(Number.isInteger(v.k) && v.k >= 2, 'k debe ser un entero ≥ 2.');
    const n = U.nMin(n => powA(v.f, v.k, n * v.k, +v.alfa), +v.potencia); return { n, N: n * v.k, real: powA(v.f, v.k, n * v.k, +v.alfa), apa: `Se necesitan ${n} participantes por grupo (N = ${n * v.k}) para detectar f = ${H.f(v.f)} con potencia de ${H.sz(H.f(+v.potencia))}.` }; };
  C['potencia-anova'] = v => { need(v.f > 0 && Number.isInteger(v.k) && v.k >= 2 && v.n >= 2, 'Revisa f, k y n.'); const p = powA(v.f, v.k, v.n * v.k, +v.alfa); return { pot: p, beta: 1 - p, lam: v.f ** 2 * v.n * v.k }; };
  const powReg = (f2, u, p, N, a) => U.powF(f2 * N, u, N - p - 1, a);
  C['tamano-muestral-regresion'] = (v, H) => { need(v.f2 > 0, 'f² debe ser mayor que 0.'); need(Number.isInteger(v.u) && v.u >= 1 && Number.isInteger(v.p) && v.p >= v.u, 'Revisa el número de predictores (los evaluados no pueden superar el total).');
    const N = U.nMin(N => powReg(v.f2, v.u, v.p, N, +v.alfa), +v.potencia, v.p + 2); return { N, real: powReg(v.f2, v.u, v.p, N, +v.alfa), apa: `Se necesitan N = ${N} participantes para detectar f² = ${H.f(v.f2)} con potencia de ${H.sz(H.f(+v.potencia))}.` }; };
  C['potencia-regresion'] = v => { need(v.f2 > 0 && v.u >= 1 && v.p >= v.u && v.n > v.p + 1, 'Revisa f², predictores y N.'); const p = powReg(v.f2, v.u, v.p, v.n, +v.alfa); return { pot: p, beta: 1 - p }; };
  const powH = (h, n, a, c) => { const zh = Math.abs(h) * Math.sqrt(n / 2), zc = U.zq(1 - a / c); return c === 2 ? U.Phi(zh - zc) + U.Phi(-zh - zc) : U.Phi(zh - zc); };
  C['tamano-muestral-proporciones'] = (v, H) => { need(v.p1 > 0 && v.p1 < 1 && v.p2 > 0 && v.p2 < 1 && v.p1 !== v.p2, 'Las proporciones deben estar entre 0 y 1 y ser distintas.');
    const h = 2 * Math.asin(Math.sqrt(v.p1)) - 2 * Math.asin(Math.sqrt(v.p2)), n = U.nMin(n => powH(h, n, +v.alfa, +v.colas), +v.potencia);
    return { h, n, N: 2 * n, apa: `Se necesitan ${n} participantes por grupo (N = ${2 * n}) para detectar ${H.f(100 * v.p1, 1)} % frente a ${H.f(100 * v.p2, 1)} %.` }; };
  C['tamano-muestral-chi'] = (v, H) => { need(v.w > 0, 'w debe ser mayor que 0.'); need(Number.isInteger(v.df) && v.df >= 1, 'gl debe ser un entero ≥ 1.');
    const N = U.nMin(N => U.powChi(v.w * v.w * N, v.df, +v.alfa), +v.potencia); return { N, real: U.powChi(v.w * v.w * N, v.df, +v.alfa), apa: `Se necesitan N = ${N} casos para detectar w = ${H.f(v.w)} con gl = ${v.df}.` }; };
  C['tamano-muestral-media-precision'] = v => { need(v.s > 0 && v.e > 0, 'La DT y el margen de error deben ser mayores que 0.'); const z = U.zq(1 - (1 - v.conf) / 2);
    let n = Math.ceil((z * v.s / v.e) ** 2); if (Number.isFinite(v.N) && v.N > 0) n = Math.ceil(n / (1 + (n - 1) / v.N)); return { n }; };
  C['tamano-muestral-proporcion-precision'] = v => { need(v.p > 0 && v.p < 1, 'p debe estar entre 0 y 1.'); need(v.e > 0 && v.e < 1, 'El margen de error debe estar entre 0 y 1 (por ejemplo 0,03).'); const z = U.zq(1 - (1 - v.conf) / 2);
    const n0 = Math.ceil(z * z * v.p * (1 - v.p) / (v.e * v.e)); return { n: n0, nf: Number.isFinite(v.N) && v.N > 0 ? Math.ceil(n0 / (1 + (n0 - 1) / v.N)) : n0 }; };
  C['efecto-minimo-detectable'] = (v, H) => { need(Number.isInteger(v.n) && v.n >= 2, 'n debe ser un entero ≥ 2.'); const tipo = v.diseno, target = +v.potencia;
    const d = U.bisect(d => U.powT(d, v.n, +v.alfa, +v.colas, tipo) - target, 1e-4, 20); return { d, apa: `Con n = ${v.n}${tipo === 'ind' ? ' por grupo' : ''}, el efecto mínimo detectable con potencia de ${H.sz(H.f(target))} es d = ${H.f(d)}.` }; };
  const powRMSEA = (N, df, e0, ea, a) => { const l0 = (N - 1) * df * e0 * e0, la = (N - 1) * df * ea * ea;
    if (ea > e0) { const cr = U.ncChiInv(1 - a, df, l0); return 1 - U.ncChiCdf(cr, df, la); } const cr = U.ncChiInv(a, df, l0); return U.ncChiCdf(cr, df, la); };
  C['potencia-sem-rmsea'] = (v, H) => { need(Number.isInteger(v.df) && v.df >= 1, 'gl del modelo debe ser un entero ≥ 1.'); need(v.e0 >= 0 && v.ea >= 0 && v.e0 !== v.ea, 'Los RMSEA deben ser no negativos y distintos.');
    const p = Number.isFinite(v.n) && v.n > 1 ? powRMSEA(v.n, v.df, v.e0, v.ea, +v.alfa) : NaN, N = U.nMin(N => powRMSEA(N, v.df, v.e0, v.ea, +v.alfa), +v.potencia, 10, 1e5);
    return { pot: p, N, apa: `Para gl = ${v.df}, ε₀ = ${H.sz(H.f(v.e0, 3))} y εₐ = ${H.sz(H.f(v.ea, 3))}, se necesitan N = ${N} casos para una potencia de ${H.sz(H.f(+v.potencia))}.` }; };

  // Metaanálisis
  C['metaanalisis-efectos'] = (v, H) => { lista(v.es, 2, 'Introduce al menos 2 tamaños del efecto.'); lista(v.vs, 2, 'Introduce las varianzas.'); need(v.es.length === v.vs.length, 'Debe haber tantas varianzas como efectos.');
    const vv = v.medida === 'ee' ? v.vs.map(s => s * s) : v.vs; need(vv.every(x => x > 0), 'Las varianzas deben ser mayores que 0.');
    const k = v.es.length, w = vv.map(x => 1 / x), sw = sum(w), fe = sum(v.es.map((y, i) => w[i] * y)) / sw, Q = sum(v.es.map((y, i) => w[i] * (y - fe) ** 2)), df = k - 1;
    const t2 = Math.max(0, (Q - df) / (sw - sum(w.map(x => x * x)) / sw)), ws = vv.map(x => 1 / (x + t2)), sws = sum(ws), re = sum(v.es.map((y, i) => ws[i] * y)) / sws, sere = Math.sqrt(1 / sws);
    const I2 = Q > df ? (Q - df) / Q : 0, pq = 1 - J.chisquare.cdf(Q, df), tpi = k > 2 ? J.studentt.inv(.975, k - 2) : NaN, pi = Math.sqrt(t2 + sere * sere);
    return { fe, re, icre: H.ci(re - 1.959964 * sere, re + 1.959964 * sere), p: U.pTwo(re / sere), Q, pq, I2: 100 * I2, t2, ip: k > 2 ? H.ci(re - tpi * pi, re + tpi * pi) : '—',
      apa: `efecto combinado = ${H.f(re)}, IC 95 % [${H.f(re - 1.959964 * sere)}; ${H.f(re + 1.959964 * sere)}], Q(${df}) = ${H.f(Q)}, ${H.apaP(pq)}, I² = ${H.f(100 * I2, 1)} %, τ² = ${H.f(t2, 3)}` }; };
  C['varianza-d'] = v => { need(fin(v.d, v.n1, v.n2) && v.n1 >= 2 && v.n2 >= 2, 'Rellena todos los campos.');
    const N = v.n1 + v.n2, vd = N / (v.n1 * v.n2) + v.d ** 2 / (2 * N), Jc = 1 - 3 / (4 * N - 9); return { vd, sed: Math.sqrt(vd), g: Jc * v.d, vg: Jc * Jc * vd, w: 1 / vd }; };
  C['stouffer'] = (v, H) => { lista(v.ps, 2, 'Introduce al menos 2 valores.'); const z = v.tipo === 'p' ? (need(v.ps.every(p => p > 0 && p < 1), 'Los p-valores deben estar entre 0 y 1.'), v.ps.map(p => U.zq(1 - p))) : v.ps;
    const k = z.length, Z = sum(z) / Math.sqrt(k), p = 1 - U.Phi(Z); return { Z, p, k, apa: `Z de Stouffer = ${H.f(Z)}, ${H.apaP(p)} (unilateral)` }; };
  C['fisher-p'] = (v, H) => { lista(v.ps, 2, 'Introduce al menos 2 p-valores.'); need(v.ps.every(p => p > 0 && p <= 1), 'Los p-valores deben estar entre 0 y 1.');
    const chi = -2 * sum(v.ps.map(Math.log)), df = 2 * v.ps.length, p = 1 - J.chisquare.cdf(chi, df); return { chi, df, p, apa: `χ²(${df}) = ${H.f(chi)}, ${H.apaP(p)}` }; };
  C['fail-safe-n'] = v => { lista(v.zs, 1, 'Introduce al menos un valor de z.'); const k = v.zs.length, za = U.zq(1 - +v.alfa), S = sum(v.zs), n = Math.max(0, S * S / (za * za) - k);
    return { nfs: Math.floor(n), umbral: 5 * k + 10, rob: n > 5 * k + 10 ? 'Robusto según el criterio 5k + 10' : 'No robusto según el criterio 5k + 10' }; };
  C['promedio-correlaciones'] = (v, H) => { lista(v.rs, 2, 'Introduce al menos 2 correlaciones.'); lista(v.ns, 2, 'Introduce los tamaños.'); need(v.rs.length === v.ns.length, 'Debe haber tantos tamaños como correlaciones.');
    need(v.rs.every(r => Math.abs(r) < 1) && v.ns.every(n => n > 3), 'Las r deben estar entre −1 y 1 y cada n ser mayor que 3.');
    const w = v.ns.map(n => n - 3), z = v.rs.map(Math.atanh), zm = sum(z.map((x, i) => w[i] * x)) / sum(w), se = 1 / Math.sqrt(sum(w)), Q = sum(z.map((x, i) => w[i] * (x - zm) ** 2));
    return { rm: Math.tanh(zm), ic: H.ci(Math.tanh(zm - 1.959964 * se), Math.tanh(zm + 1.959964 * se)), Q, pq: 1 - J.chisquare.cdf(Q, v.rs.length - 1), rsimple: mean(v.rs) }; };

  // Bayesiana
  C['factor-bayes-t'] = (v, H) => { need(fin(v.t, v.n1) && v.n1 >= 2, 'Rellena t y n.'); need(v.rscale > 0, 'La escala del prior debe ser mayor que 0.');
    let N, nu; if (v.diseno === 'ind') { need(Number.isFinite(v.n2) && v.n2 >= 2, 'Introduce n₂.'); N = v.n1 * v.n2 / (v.n1 + v.n2); nu = v.n1 + v.n2 - 2; } else { N = v.n1; nu = v.n1 - 1; }
    const r2 = v.rscale ** 2, lik0 = -(nu + 1) / 2 * Math.log(1 + v.t * v.t / nu);
    const f = g => -0.5 * Math.log(1 + N * g) - (nu + 1) / 2 * Math.log(1 + v.t * v.t / ((1 + N * g) * nu)) + 0.5 * Math.log(r2 / (2 * Math.PI)) - 1.5 * Math.log(g) - r2 / (2 * g);
    const a = -25, b = 25, m = 4000, hh = (b - a) / m; let s = 0;
    for (let i = 0; i <= m; i++) { const u = a + i * hh, g = Math.exp(u), w = i === 0 || i === m ? 1 : (i % 2 ? 4 : 2); s += w * Math.exp(f(g) - lik0 + u); }
    const bf = s * hh / 3; return { bf10: bf, bf01: 1 / bf, lbf: Math.log(bf), post: bf / (1 + bf), apa: `BF₁₀ = ${H.f(bf)}` }; };
  C['factor-bayes-bic'] = (v, H) => { need(fin(v.bic0, v.bic1), 'Introduce los dos BIC.'); const bf = Math.exp((v.bic0 - v.bic1) / 2); return { bf10: bf, bf01: 1 / bf, post: bf / (1 + bf), apa: `BF₁₀ ≈ ${H.f(bf)}` }; };
  C['beta-binomial'] = (v, H) => { need(v.a > 0 && v.b > 0, 'Los parámetros del prior deben ser mayores que 0.'); need(Number.isInteger(v.x) && Number.isInteger(v.n) && v.x >= 0 && v.x <= v.n, 'x y n deben ser enteros con 0 ≤ x ≤ n.');
    const A = v.a + v.x, B = v.b + v.n - v.x, pc = Math.round(v.conf * 100), al = 1 - v.conf, lo = J.beta.inv(al / 2, A, B), hi = J.beta.inv(1 - al / 2, A, B), th = Number.isFinite(v.ref) ? v.ref : 0.5;
    return { post: `Beta(${H.f(A, 1)}; ${H.f(B, 1)})`, media: A / (A + B), moda: A > 1 && B > 1 ? (A - 1) / (A + B - 2) : NaN, ic: H.ci(lo, hi, 3), pmay: 1 - J.beta.cdf(th, A, B), apa: `media posterior = ${H.f(A / (A + B), 3)}, intervalo de credibilidad ${pc} % [${H.f(lo, 3)}; ${H.f(hi, 3)}]` }; };
  C['normal-normal'] = (v, H) => { need(fin(v.m0, v.s0, v.xm, v.sigma, v.n) && v.s0 > 0 && v.sigma > 0 && v.n >= 1, 'Rellena todos los campos (DT > 0).');
    const p0 = 1 / v.s0 ** 2, pd = v.n / v.sigma ** 2, pp = p0 + pd, mp = (p0 * v.m0 + pd * v.xm) / pp, sp = Math.sqrt(1 / pp), z = U.zq(1 - (1 - v.conf) / 2);
    return { mp, sp, ic: H.ci(mp - z * sp, mp + z * sp), peso: pd / pp }; };
  C['probabilidad-posterior'] = v => { need(v.bf > 0, 'El factor de Bayes debe ser mayor que 0.'); need(v.prior > 0 && v.prior < 1, 'La probabilidad previa debe estar entre 0 y 1.');
    const po = v.bf * v.prior / (1 - v.prior); return { odds: po, post: po / (1 + po), prioro: v.prior / (1 - v.prior) }; };
})(typeof window !== 'undefined' ? window : globalThis);

/* Interfaz común: lee los campos de la ficha, llama a CALC[slug] y muestra los resultados. */
(function () {
  if (typeof document === 'undefined') return;
  const $ = s => document.querySelector(s);
  const f = (x, k = 2) => { if (!Number.isFinite(x)) return x === Infinity ? '∞' : '—'; const s = x.toLocaleString('es-ES', { minimumFractionDigits: k, maximumFractionDigits: k }); return /^-0(,0+)?$/.test(s) ? s.slice(1) : s; };
  const sz = s => String(s).replace(/^(-?)0,/, '$1,');
  const H = { f, sz, apaP: p => p < 0.001 ? 'p < ,001' : 'p = ' + sz(f(p, 3)), ci: (a, b, k = 2) => `[${f(a, k)}; ${f(b, k)}]` };
  const auto = x => { if (!Number.isFinite(x)) return x === Infinity ? '∞' : '—'; const a = Math.abs(x);
    if (Number.isInteger(x) && a < 1e15) return x.toLocaleString('es-ES');
    if (a >= 1e9 || (a < 1e-4 && x !== 0)) return x.toExponential(4).replace('.', ',');
    return x.toLocaleString('es-ES', { maximumFractionDigits: 6 }); };
  const fmt = (x, t) => { if (x === undefined || x === null) return '—'; if (typeof x === 'string') return x;
    if (t === 'p') return !Number.isFinite(x) ? '—' : x < 0.001 ? '< 0,001' : f(x, 3);
    if (t === 'i') return Number.isFinite(x) ? Math.round(x).toLocaleString('es-ES') : '—';
    if (t === 'a') return auto(x); if (t === 't') return Number.isFinite(x) ? f(x, 2) : '—';
    return f(x, +t); };
  const num = s => { s = s.trim(); if (s === '') return NaN; return Number(s.replace(/\s/g, '').replace(',', '.')); };
  const toks = s => s.split(/[\s;\t]+/).filter(Boolean).map(t => Number(t.replace(',', '.')));
  const leer = el => { const t = el.dataset.t, val = el.value;
    if (t === 'n') return num(val); if (t === 's') return val;
    if (t === 'l') return val.trim() === '' ? [] : toks(val);
    return val.trim().split(/\n+/).filter(l => l.trim()).map(toks); };

  function iniciar() {
    const main = document.querySelector('main[data-calc]'); if (!main) return;
    const fn = window.CALC && window.CALC[main.dataset.calc]; if (!fn) return;
    const card = main.querySelector('.card'), err = $('#err'), apaBox = $('#apa-box'), outs = [...main.querySelectorAll('[data-out]')];
    const run = () => {
      const v = {}; card.querySelectorAll('[data-k]').forEach(el => { v[el.dataset.k] = leer(el); });
      try {
        const r = fn(v, H); err.hidden = true;
        outs.forEach(o => { o.textContent = fmt(r[o.dataset.out], o.dataset.f); });
        if (r.apa) { $('#apa').textContent = r.apa; apaBox.hidden = false; } else apaBox.hidden = true;
      } catch (e) { err.textContent = e.message; err.hidden = false; outs.forEach(o => o.textContent = '—'); apaBox.hidden = true; }
    };
    card.addEventListener('input', run); card.addEventListener('change', run); run();
  }

  document.addEventListener('click', e => {
    const b = e.target.closest('[data-copy]'); if (!b) return;
    const el = document.getElementById(b.dataset.copy), txt = el.textContent;
    const ok = () => { b.textContent = 'Copiado'; setTimeout(() => b.textContent = 'Copiar', 1500); };
    const sel = () => { const r = document.createRange(); r.selectNodeContents(el); const s = getSelection(); s.removeAllRanges(); s.addRange(r); };
    if (navigator.clipboard) navigator.clipboard.writeText(txt).then(ok, sel); else sel();
  });

  // Catálogo: pestañas de categoría y buscador
  function catalogo() {
    const tabs = document.querySelectorAll('.cat-tabs [data-cat]'), q = document.getElementById('buscar'); if (!tabs.length) return;
    const norm = s => s.toLowerCase().normalize('NFD').replace(/[̀-ͯ]/g, '');
    let actual = 'todas';
    const aplicar = () => { const t = q ? norm(q.value.trim()) : '';
      document.querySelectorAll('.cat').forEach(sec => { let vis = 0;
        sec.querySelectorAll('.ccard').forEach(c => { const ok = (actual === 'todas' || sec.dataset.cat === actual) && (!t || norm(c.textContent).includes(t)); c.hidden = !ok; if (ok) vis++; });
        sec.hidden = vis === 0; });
      tabs.forEach(b => b.setAttribute('aria-selected', b.dataset.cat === actual ? 'true' : 'false')); };
    tabs.forEach(b => b.addEventListener('click', () => { actual = b.dataset.cat; if (history.replaceState) history.replaceState(null, '', actual === 'todas' ? location.pathname : '#' + actual); aplicar(); }));
    if (q) q.addEventListener('input', aplicar);
    const h = location.hash.slice(1); if (h && document.querySelector(`.cat-tabs [data-cat="${h}"]`)) actual = h;
    aplicar();
  }

  const start = () => { iniciar(); catalogo(); const cur = document.querySelector('.cat-nav [aria-current]'); if (cur) cur.parentElement.scrollLeft = cur.offsetLeft - 24; };
  if (document.readyState === 'loading') document.addEventListener('DOMContentLoaded', start); else start();
})();

document.querySelectorAll('[data-url]').forEach(e => e.textContent = location.href.split('#')[0]);
