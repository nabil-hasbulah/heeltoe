/* Heel & Toe analysis engine. Pure functions, no DOM. */
var HT = (function () {
  var G = 9.80665, DT = 0.1;

  function clamp(x, a, b) { return x < a ? a : x > b ? b : x; }
  function score(x, good, bad) { return Math.round(clamp(100 * (bad - x) / (bad - good), 0, 100)); }
  function mean(a) { var s = 0; for (var i = 0; i < a.length; i++) s += a[i]; return a.length ? s / a.length : 0; }
  function median(a) { var b = a.slice().sort(function (x, y) { return x - y; }); return b.length ? b[b.length >> 1] : 0; }

  // centred box filter, applied twice (triangular), window in samples
  function smooth(a, w) {
    w = Math.max(1, Math.round(w));
    if (w < 2) return a.slice();
    var out = box(box(a, w), w);
    return out;
  }
  function box(a, w) {
    var n = a.length, out = new Array(n), h = w >> 1, s = 0, c = 0, lo = 0, hi = -1;
    for (var i = 0; i < n; i++) {
      var L = Math.max(0, i - h), R = Math.min(n - 1, i + h);
      while (hi < R) { hi++; s += a[hi]; c++; }
      while (lo < L) { s -= a[lo]; lo++; c--; }
      out[i] = s / c;
    }
    return out;
  }
  function deriv(a, dt) {
    var n = a.length, out = new Array(n);
    for (var i = 0; i < n; i++) {
      var p = Math.max(0, i - 1), q = Math.min(n - 1, i + 1);
      out[i] = q > p ? (a[q] - a[p]) / ((q - p) * dt) : 0;
    }
    return out;
  }

  /* ---------- parsers: both return {startedAt, gps:[[ms,lat,lon,spd,hdg,acc,alt]], imu:[[ms,ax,ay,az]]|null, source} ---------- */
  function parseGPX(text) {
    var pts = [], re = /<trkpt\b([^>]*?)(\/>|>([\s\S]*?)<\/trkpt>)/g, m;
    while ((m = re.exec(text))) {
      var a = m[1], body = m[3] || '';
      var lat = parseFloat((/lat\s*=\s*["']([^"']+)/.exec(a) || [])[1]);
      var lon = parseFloat((/lon\s*=\s*["']([^"']+)/.exec(a) || [])[1]);
      var tm = (/<time>([^<]+)<\/time>/.exec(body) || [])[1];
      var sp = (/<(?:[\w]+:)?speed>([^<]+)</.exec(body) || [])[1];
      var el = (/<ele>([^<]+)<\/ele>/.exec(body) || [])[1];
      var t = tm ? Date.parse(tm) : NaN;
      if (isFinite(lat) && isFinite(lon) && isFinite(t)) pts.push([t, lat, lon, sp != null && isFinite(+sp) ? +sp : null, null, null, el != null ? +el : null]);
    }
    if (pts.length < 2) throw new Error('No timed track points found in this GPX file. Export it with timestamps on.');
    pts.sort(function (x, y) { return x[0] - y[0]; });
    var t0 = pts[0][0];
    return { startedAt: t0, gps: pts.map(function (p) { return [p[0] - t0].concat(p.slice(1)); }), imu: null, source: 'gpx' };
  }
  function parseRecording(obj) {
    if (!obj || obj.format !== 'heeltoe-rec' || !Array.isArray(obj.gps)) throw new Error('This file is not a Heel & Toe recording.');
    return { startedAt: obj.startedAt || Date.now(), gps: obj.gps, imu: Array.isArray(obj.imu) && obj.imu.length > 50 ? obj.imu : null, source: 'recorder', device: obj.device || '' };
  }
  function parseFile(name, text) {
    var t = text.trim();
    if (t[0] === '{') return parseRecording(JSON.parse(t));
    if (/<gpx[\s>]/i.test(t)) return parseGPX(t);
    throw new Error('Unsupported file. Use a Heel & Toe recording (.json) or a .gpx track.');
  }

  /* ---------- core processing ---------- */
  function process(rec) {
    var raw = rec.gps.filter(function (g) { return isFinite(g[1]) && isFinite(g[2]) && (g[5] == null || g[5] <= 30); });
    raw.sort(function (a, b) { return a[0] - b[0]; });
    var gps = [];
    for (var i = 0; i < raw.length; i++) if (!gps.length || raw[i][0] - gps[gps.length - 1][0] >= 50) gps.push(raw[i]);
    if (gps.length < 20) throw new Error('Too few GPS points (' + gps.length + '). Record at least a minute of driving with a clear sky view.');

    var lat0 = gps[0][1], lon0 = gps[0][2], kx = 111320 * Math.cos(lat0 * Math.PI / 180), ky = 110540;
    var gt = gps.map(function (g) { return g[0] / 1000; });
    var intervals = []; for (i = 1; i < gt.length; i++) intervals.push(gt[i] - gt[i - 1]);
    var gpsInt = median(intervals);
    var T0 = gt[0], N = Math.floor((gt[gt.length - 1] - T0) / DT) + 1;
    if (N > 10 * 3600 * 4) throw new Error('Recording is longer than 4 hours. Split it into shorter drives.');

    var hasSpd = gps.filter(function (g) { return g[3] != null && g[3] >= 0; }).length > gps.length * 0.8;
    var X = new Array(N), Y = new Array(N), S = new Array(N), Al = new Array(N);
    var j = 0;
    for (i = 0; i < N; i++) {
      var t = T0 + i * DT;
      while (j < gt.length - 2 && gt[j + 1] < t) j++;
      var a = gps[j], b = gps[j + 1], f = clamp((t - gt[j]) / Math.max(1e-6, gt[j + 1] - gt[j]), 0, 1);
      X[i] = ((a[2] + (b[2] - a[2]) * f) - lon0) * kx;
      Y[i] = ((a[1] + (b[1] - a[1]) * f) - lat0) * ky;
      var sa = a[3] != null && a[3] >= 0 ? a[3] : null, sb = b[3] != null && b[3] >= 0 ? b[3] : null;
      S[i] = sa != null && sb != null ? sa + (sb - sa) * f : (sa != null ? sa : sb);
      Al[i] = a[6] != null && b[6] != null ? a[6] + (b[6] - a[6]) * f : null;
    }
    // smoothing window adapts to GPS rate
    var pw = Math.max(5, Math.round(Math.max(0.6, gpsInt * 1.2) / DT));
    var Xs = smooth(X, pw), Ys = smooth(Y, pw);
    var vx = deriv(Xs, DT), vy = deriv(Ys, DT);
    var vPos = vx.map(function (x, k) { return Math.hypot(x, vy[k]); });
    var V;
    if (hasSpd) { V = S.map(function (s, k) { return s == null ? vPos[k] : s; }); V = smooth(V, 4); }
    else V = vPos;

    var hdg = vx.map(function (x, k) { return Math.atan2(vy[k], x); });
    for (i = 1; i < N; i++) { var d = hdg[i] - hdg[i - 1]; while (d > Math.PI) { hdg[i] -= 2 * Math.PI; d -= 2 * Math.PI; } while (d < -Math.PI) { hdg[i] += 2 * Math.PI; d += 2 * Math.PI; } }
    var yaw = deriv(hdg, DT).map(function (w, k) { return vPos[k] > 2 ? w : 0; });
    yaw = smooth(yaw, 6);
    var latG = V.map(function (v, k) { return clamp(v * yaw[k] / G, -1.5, 1.5); });
    var lonG = smooth(deriv(V, DT), hasSpd ? 6 : pw).map(function (x) { return clamp(x / G, -1.5, 1.5); });

    // ---------- IMU alignment + complementary fusion ----------
    var imuInfo = { used: false };
    var fLat = latG, fLon = lonG;
    if (rec.imu) {
      var r = alignImu(rec.imu, T0 * 1000, N, latG, lonG, V);
      imuInfo = r.info;
      if (r.ok) {
        var LP = 25; // 2.5 s
        var gl = smooth(latG, LP), go = smooth(lonG, LP), il = smooth(r.lat, LP), io = smooth(r.lon, LP);
        fLat = r.lat.map(function (x, k) { return gl[k] + x - il[k]; });
        fLon = r.lon.map(function (x, k) { return go[k] + x - io[k]; });
        fLat = smooth(fLat, 2); fLon = smooth(fLon, 2);
      }
    }
    var jw = 3; // light smoothing before jerk
    var sl = smooth(fLat, jw), so = smooth(fLon, jw);
    var dLat = deriv(sl, DT), dLon = deriv(so, DT);
    var jerk = dLat.map(function (x, k) { return Math.hypot(x, dLon[k]); });
    var moving = V.map(function (v) { return v > 3; });

    // ---------- trim stationary ends ----------
    var first = 0, last = N - 1;
    while (first < N && !moving[first]) first++;
    while (last > first && !moving[last]) last--;
    first = Math.max(0, first - 20); last = Math.min(N - 1, last + 20);
    if (last - first < 300) throw new Error('Less than 30 seconds of movement in this recording.');
    function cut(a) { return a.slice(first, last + 1); }
    var D = { x: cut(Xs), y: cut(Ys), v: cut(V), lat: cut(fLat), lon: cut(fLon), jerk: cut(jerk), dLat: cut(dLat), dLon: cut(dLon), yaw: cut(yaw), alt: cut(Al), moving: cut(moving) };
    var n = D.v.length;

    // ---------- summary ----------
    var dist = 0, movT = 0, maxV = 0, js = [];
    for (i = 1; i < n; i++) dist += Math.hypot(D.x[i] - D.x[i - 1], D.y[i] - D.y[i - 1]);
    for (i = 0; i < n; i++) if (D.moving[i]) { movT += DT; js.push(D.jerk[i] * D.jerk[i]); if (D.v[i] > maxV) maxV = D.v[i]; }
    var rmsJerk = Math.sqrt(mean(js));

    var corners = findCorners(D);
    var scores = overallScores(corners, rmsJerk);

    // time-ordered elevation gain (only if altitude present)
    var climb = null;
    if (D.alt.filter(function (a) { return a != null; }).length > n * 0.8) {
      var as = smooth(D.alt.map(function (a, k) { return a == null ? (D.alt[k - 1] || 0) : a; }), 50); climb = 0;
      for (i = 10; i < n; i += 10) if (as[i] > as[i - 10]) climb += as[i] - as[i - 10];
      climb = Math.round(climb);
    }

    return {
      origin: [lat0, lon0], startedAt: rec.startedAt + (T0 + first * DT) * 1000, dt: DT, n: n, D: D,
      summary: { dist: Math.round(dist), moveSec: Math.round(movT), totalSec: Math.round(n * DT), maxKmh: Math.round(maxV * 3.6), avgKmh: movT ? Math.round(dist / movT * 3.6) : 0, rmsJerk: +rmsJerk.toFixed(3), climb: climb, gpsHz: +(1 / gpsInt).toFixed(1) },
      imu: imuInfo, corners: corners, scores: scores, source: rec.source
    };
  }

  function alignImu(imu, t0ms, N, latG, lonG, V) {
    // bin to the 10 Hz grid
    var sx = new Array(N).fill(0), sy = new Array(N).fill(0), sz = new Array(N).fill(0), c = new Array(N).fill(0);
    for (var i = 0; i < imu.length; i++) {
      var k = Math.round((imu[i][0] - t0ms) / 100);
      if (k < 0 || k >= N) continue;
      sx[k] += imu[i][1]; sy[k] += imu[i][2]; sz[k] += imu[i][3]; c[k]++;
    }
    var filled = c.filter(function (x) { return x > 0; }).length;
    if (filled < N * 0.5) return { ok: false, info: { used: false, reason: 'Motion data covers too little of the drive' } };
    var ax = [], ay = [], az = [], last = null;
    for (k = 0; k < N; k++) {
      if (c[k]) last = [sx[k] / c[k], sy[k] / c[k], sz[k] / c[k]];
      var v = last || [0, 0, 0]; ax.push(v[0]); ay.push(v[1]); az.push(v[2]);
    }
    var m = [mean(ax), mean(ay), mean(az)], mg = Math.hypot(m[0], m[1], m[2]);
    if (mg < 5) return { ok: false, info: { used: false, reason: 'Motion data has no gravity reference' } };
    var gh = [m[0] / mg, m[1] / mg, m[2] / mg];
    var ref = Math.abs(gh[0]) < 0.8 ? [1, 0, 0] : [0, 1, 0];
    var e1 = norm(cross(gh, ref)), e2 = cross(gh, e1);
    var P = [], Q = [];
    for (k = 0; k < N; k++) {
      var l = [ax[k] - m[0], ay[k] - m[1], az[k] - m[2]];
      P.push(dot(l, e1) / G); Q.push(dot(l, e2) / G);
    }
    P = smooth(P, 3); Q = smooth(Q, 3);
    var Ps = smooth(P, 6), Qs = smooth(Q, 6); // match the GPS bandwidth for fitting
    var best = null;
    for (var L = -25; L <= 25; L++) {
      var fo = fit(lonG, Ps, Qs, V, L), fl = fit(latG, Ps, Qs, V, L);
      if (!fo || !fl) continue;
      var s = fo.r2 + fl.r2;
      if (!best || s > best.s) best = { s: s, L: L, fo: fo, fl: fl };
    }
    if (!best) return { ok: false, info: { used: false, reason: 'Not enough moving data to align the phone' } };
    var fwd = norm2([best.fo.a, best.fo.b]);
    var lat1 = [-fwd[1], fwd[0]];
    if (lat1[0] * best.fl.a + lat1[1] * best.fl.b < 0) lat1 = [-lat1[0], -lat1[1]];
    var info = { used: false, lagSec: +(best.L * DT).toFixed(1), fitLon: +best.fo.r2.toFixed(2), fitLat: +best.fl.r2.toFixed(2) };
    if (best.fo.r2 < 0.2 || best.fl.r2 < 0.35) { info.reason = 'Phone movement did not match the GPS track (loose mount?)'; return { ok: false, info: info }; }
    var oLat = new Array(N), oLon = new Array(N);
    for (k = 0; k < N; k++) {
      var q = clamp(k + best.L, 0, N - 1);
      oLon[k] = clamp(fwd[0] * P[q] + fwd[1] * Q[q], -1.5, 1.5);
      oLat[k] = clamp(lat1[0] * P[q] + lat1[1] * Q[q], -1.5, 1.5);
    }
    info.used = true;
    return { ok: true, lat: oLat, lon: oLon, info: info };
  }
  function fit(y, P, Q, V, L) {
    var spp = 0, sqq = 0, spq = 0, syp = 0, syq = 0, syy = 0, n = 0, N = y.length;
    for (var k = 0; k < N; k += 2) {
      var q = k + L; if (q < 0 || q >= N || V[k] < 4) continue;
      var p = P[q], r = Q[q], yy = y[k];
      spp += p * p; sqq += r * r; spq += p * r; syp += yy * p; syq += yy * r; syy += yy * yy; n++;
    }
    if (n < 200 || syy <= 0) return null;
    var det = spp * sqq - spq * spq; if (Math.abs(det) < 1e-12) return null;
    var a = (syp * sqq - syq * spq) / det, b = (syq * spp - syp * spq) / det;
    var sse = syy - 2 * (a * syp + b * syq) + a * a * spp + 2 * a * b * spq + b * b * sqq;
    return { a: a, b: b, r2: Math.max(0, 1 - sse / syy) };
  }
  function cross(a, b) { return [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]]; }
  function dot(a, b) { return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]; }
  function norm(a) { var m = Math.hypot(a[0], a[1], a[2]) || 1; return [a[0] / m, a[1] / m, a[2] / m]; }
  function norm2(a) { var m = Math.hypot(a[0], a[1]) || 1; return [a[0] / m, a[1] / m]; }

  /* ---------- corners ---------- */
  var LAT_ON = 0.18, BRAKE = -0.08, THR = 0.05;
  function findCorners(D) {
    var n = D.v.length, runs = [], cur = null;
    for (var i = 0; i < n; i++) {
      var s = Math.abs(D.lat[i]) > LAT_ON && D.v[i] > 4 ? Math.sign(D.lat[i]) : 0;
      if (s && cur && cur.sg === s && i - cur.e <= 8) cur.e = i;
      else if (s) { if (cur) runs.push(cur); cur = { s: i, e: i, sg: s }; }
    }
    if (cur) runs.push(cur);
    var out = [], prevEnd = 0;
    runs.forEach(function (r) {
      var dur = (r.e - r.s) * DT;
      var turn = 0; for (var k = r.s; k <= r.e; k++) turn += D.yaw[k] * DT;
      turn = Math.abs(turn) * 180 / Math.PI;
      if (dur < 1.2 || turn < 20) return;
      var c = analyseCorner(D, r.s, r.e, r.sg, prevEnd);
      c.turn = Math.round(turn);
      out.push(c); prevEnd = r.e;
    });
    out.forEach(function (c, i) { c.no = i + 1; });
    return out;
  }
  function analyseCorner(D, s, e, sg, prevEnd) {
    var n = D.v.length, apex = s, k;
    for (k = s; k <= e; k++) if (D.v[k] < D.v[apex]) apex = k;
    var peakLat = 0; for (k = s; k <= e; k++) peakLat = Math.max(peakLat, Math.abs(D.lat[k]));
    // braking block that ends closest before the apex
    var from = Math.max(prevEnd, s - 50), bEnd = -1, bStart = -1;
    for (k = apex; k >= from; k--) if (D.lon[k] < BRAKE) { bEnd = k; break; }
    if (bEnd >= 0) { bStart = bEnd; var gap = 0; for (k = bEnd; k >= from; k--) { if (D.lon[k] < BRAKE) { bStart = k; gap = 0; } else if (++gap > 3) break; } }
    var c = { s: s, e: e, apex: apex, dir: sg > 0 ? 'L' : 'R', dur: +((e - s) * DT).toFixed(1), entryKmh: Math.round(D.v[s] * 3.6), minKmh: Math.round(D.v[apex] * 3.6), exitKmh: Math.round(D.v[e] * 3.6), peakLat: +peakLat.toFixed(2) };
    c.radius = peakLat > 0.05 ? Math.round(D.v[apex] * D.v[apex] / (peakLat * G)) : null;
    var subs = [];
    if (bStart >= 0) {
      var pb = 0, on = 0, rel = 0, trail = 0;
      for (k = bStart; k <= bEnd; k++) pb = Math.min(pb, D.lon[k]);
      for (k = Math.max(0, bStart - 5); k <= Math.min(n - 1, bStart + 10); k++) on = Math.max(on, -D.dLon[k]);
      for (k = Math.max(0, bEnd - 8); k <= Math.min(n - 1, bEnd + 3); k++) rel = Math.max(rel, D.dLon[k]);
      for (k = Math.max(bStart, s); k <= bEnd; k++) if (Math.abs(D.lat[k]) > LAT_ON) trail += DT;
      c.brake = { peak: +(-pb).toFixed(2), onset: +on.toFixed(2), release: +rel.toFixed(2), trail: +trail.toFixed(1), start: bStart, end: bEnd };
      c.sBrake = score(on, 0.6, 2.5); c.sRelease = score(rel, 0.6, 2.5);
      subs.push(c.sBrake, c.sRelease);
    }
    var ti = 0; for (k = Math.max(0, s - 5); k <= Math.min(n - 1, s + 15); k++) ti = Math.max(ti, Math.abs(D.dLat[k]));
    c.turnIn = +ti.toFixed(2); c.sTurnIn = score(ti, 0.4, 1.6); subs.push(c.sTurnIn);
    var thr = -1; for (k = apex; k <= Math.min(n - 1, e + 30); k++) if (D.lon[k] > THR) { thr = k; break; }
    c.throttleAt = thr;
    if (thr >= 0) { var tj = 0; for (k = thr; k <= Math.min(n - 1, thr + 15); k++) tj = Math.max(tj, D.dLon[k]); c.throttle = +tj.toFixed(2); c.sThrottle = score(tj, 0.35, 1.4); subs.push(c.sThrottle); }
    var cs = Math.max(s, bEnd + 1), ce = thr >= 0 ? Math.min(thr, e) : e, coast = 0;
    for (k = cs; k <= ce; k++) if (D.lon[k] < -0.015 && D.lon[k] >= BRAKE) coast += DT;
    c.coast = +coast.toFixed(1); c.coastFrac = +(coast / Math.max(0.5, c.dur)).toFixed(2); c.sFlow = score(c.coastFrac, 0.08, 0.4); subs.push(c.sFlow);
    var js = []; for (k = Math.max(0, s - 10); k <= Math.min(n - 1, e + 5); k++) js.push(D.jerk[k] * D.jerk[k]);
    c.rmsJerk = +Math.sqrt(mean(js)).toFixed(2); c.sJerk = score(c.rmsJerk, 0.2, 0.8); subs.push(c.sJerk);
    c.score = Math.round(mean(subs));
    return c;
  }
  function avgOf(cs, key) { var v = cs.filter(function (c) { return c[key] != null; }).map(function (c) { return c[key]; }); return v.length ? Math.round(mean(v)) : null; }
  function overallScores(cs, rmsJerk) {
    var o = {
      braking: (function () { var v = []; cs.forEach(function (c) { if (c.brake) v.push(c.sBrake, c.sRelease); }); return v.length ? Math.round(mean(v)) : null; })(),
      turnIn: avgOf(cs, 'sTurnIn'), throttle: avgOf(cs, 'sThrottle'), flow: avgOf(cs, 'sFlow'), steady: score(rmsJerk, 0.15, 0.6)
    };
    var w = { braking: 1, turnIn: 1, throttle: 1, flow: 0.7, steady: 1 }, s = 0, ws = 0;
    for (var k in w) if (o[k] != null) { s += o[k] * w[k]; ws += w[k]; }
    o.overall = ws ? Math.round(s / ws) : null;
    return o;
  }

  /* ---------- coaching notes (rule based) ---------- */
  var SUB = {
    braking: ['Braking', 'How progressively you squeeze on and ease off the brakes.'],
    turnIn: ['Turn-in', 'How gradually cornering load builds as you steer in.'],
    throttle: ['Throttle', 'How gently you feed power back in after the slowest point.'],
    flow: ['Flow', 'How little time you spend coasting mid-corner, neither braking nor driving.'],
    steady: ['Composure', 'How calm the car is overall: fewer sudden changes in load.']
  };
  function insights(res, prev) {
    var out = [], s = res.scores, cs = res.corners;
    if (!cs.length) { out.push({ k: 'info', t: 'No corners detected', b: 'This drive had no turns sharp enough to analyse. Twisty roads give the most useful feedback.' }); return out; }
    var keys = ['braking', 'turnIn', 'throttle', 'flow', 'steady'].filter(function (k) { return s[k] != null; });
    keys.sort(function (a, b) { return s[a] - s[b]; });
    var w = keys[0];
    var tips = {
      braking: 'Squeeze the brake pedal over about half a second instead of stabbing it, and bleed the pressure off as you start to steer. A smooth release keeps weight on the front tyres when they need it.',
      turnIn: 'Turn the wheel in one steady motion, a little slower than feels natural. If you need a second correction mid-corner, your entry speed or turn-in point was off.',
      throttle: 'Pick up the throttle gently at the slowest point and add more only as you unwind the wheel. Hands and right foot should move together, in opposite directions.',
      flow: 'You are coasting through a lot of corners. Either carry a little brake into the turn (trail braking) or get to a light maintenance throttle sooner, so the car always has a settled load on it.',
      steady: 'The car is changing load often. Look further ahead so your inputs can be earlier and smaller.'
    };
    out.push({ k: 'focus', t: 'Work on: ' + SUB[w][0].toLowerCase(), b: tips[w] });
    var best = keys[keys.length - 1];
    if (best !== w && s[best] >= 70) out.push({ k: 'good', t: 'Strongest area: ' + SUB[best][0].toLowerCase(), b: 'Scored ' + s[best] + ' this drive. Keep it as your reference feeling.' });
    var L = cs.filter(function (c) { return c.dir === 'L'; }), R = cs.filter(function (c) { return c.dir === 'R'; });
    if (L.length >= 3 && R.length >= 3) {
      var sl = mean(L.map(function (c) { return c.score; })), sr = mean(R.map(function (c) { return c.score; }));
      if (Math.abs(sl - sr) >= 8) out.push({ k: 'pattern', t: (sl < sr ? 'Left' : 'Right') + '-handers are rougher', b: 'Average corner score ' + Math.round(Math.min(sl, sr)) + ' vs ' + Math.round(Math.max(sl, sr)) + '. Look at your turn-in and braking on that side specifically.' });
    }
    var braked = cs.filter(function (c) { return c.brake; });
    if (braked.length >= 4) {
      var trail = braked.filter(function (c) { return c.brake.trail >= 0.4; }).length / braked.length;
      if (trail < 0.25) out.push({ k: 'pattern', t: 'Brakes off before turning in', b: 'Only ' + Math.round(trail * 100) + '% of braked corners carry some brake into the turn. Trailing a little brake pressure helps the front bite on entry.' });
      var abrupt = braked.filter(function (c) { return c.sRelease < 45; }).length / braked.length;
      if (abrupt > 0.4) out.push({ k: 'pattern', t: 'Sudden brake releases', b: Math.round(abrupt * 100) + '% of brake releases were abrupt. Letting go too fast shifts weight to the rear just as you turn in.' });
    }
    if (prev && prev.scores && prev.scores.overall != null && s.overall != null) {
      var d = s.overall - prev.scores.overall;
      out.push({ k: d >= 0 ? 'good' : 'pattern', t: d >= 0 ? 'Smoother than last time' : 'Rougher than last time', b: 'Overall ' + s.overall + ' vs ' + prev.scores.overall + ' on your previous drive on this road.' });
    }
    return out;
  }

  /* ---------- storage packing (5 Hz) ---------- */
  function pack(res) {
    var D = res.D, step = 2, cols = { x: [], y: [], v: [], lat: [], lon: [], j: [] };
    for (var i = 0; i < res.n; i += step) {
      cols.x.push(Math.round(D.x[i] * 10)); cols.y.push(Math.round(D.y[i] * 10));
      cols.v.push(Math.round(D.v[i] * 36)); cols.lat.push(Math.round(D.lat[i] * 100));
      cols.lon.push(Math.round(D.lon[i] * 100)); cols.j.push(Math.round(D.jerk[i] * 100));
    }
    return cols;
  }
  function unpack(cols) {
    return { x: cols.x.map(function (v) { return v / 10; }), y: cols.y.map(function (v) { return v / 10; }), v: cols.v.map(function (v) { return v / 36; }), lat: cols.lat.map(function (v) { return v / 100; }), lon: cols.lon.map(function (v) { return v / 100; }), jerk: cols.j.map(function (v) { return v / 100; }) };
  }
  // compact corner record (indices at 5 Hz)
  function slimCorner(c) {
    var o = {}; for (var k in c) o[k] = c[k];
    ['s', 'e', 'apex', 'throttleAt'].forEach(function (k) { if (o[k] >= 0) o[k] = Math.round(o[k] / 2); });
    if (o.brake) { o.brake = Object.assign({}, o.brake, { start: Math.round(o.brake.start / 2), end: Math.round(o.brake.end / 2) }); }
    return o;
  }

  /* ---------- synthetic example drive ---------- */
  function rng(seed) { return function () { seed = (seed * 1664525 + 1013904223) >>> 0; return seed / 4294967296; }; }
  function gauss(r) { return Math.sqrt(-2 * Math.log(r() + 1e-9)) * Math.cos(2 * Math.PI * r()); }
  function synth(opts) {
    opts = opts || {};
    var r = rng(opts.seed || 7);
    // road: [length m, radius m (signed, + left, 0 straight)]
    var road = opts.road || [[250, 0], [110, 70], [120, 0], [90, -45], [60, 0], [140, 120], [200, 0], [80, -35], [40, 0], [80, 40], [300, 0], [160, -90], [90, 0], [120, 55], [70, -60], [260, 0], [100, 30], [150, 0], [130, -80], [220, 0], [90, 50], [60, 0], [110, -40], [350, 0], [140, 100], [80, 0], [95, -50], [200, 0]];
    var ds = 1, kap = [];
    road.forEach(function (sgm) { for (var i = 0; i < sgm[0]; i++) kap.push(sgm[1] ? 1 / sgm[1] : 0); });
    kap = smooth(kap, 25);
    var M = kap.length, hd = 0, px = [0], py = [0], hs = [0];
    for (var i = 1; i < M; i++) { hd += kap[i] * ds; hs.push(hd); px.push(px[i - 1] + Math.cos(hd) * ds); py.push(py[i - 1] + Math.sin(hd) * ds); }
    // speed profile; some corners driven rough
    var vmax = kap.map(function (k) { return Math.min(27, Math.sqrt(0.5 * G / Math.max(Math.abs(k), 1e-4))); });
    var v = vmax.slice(); v[0] = 3;
    for (i = 1; i < M; i++) v[i] = Math.min(v[i], Math.sqrt(v[i - 1] * v[i - 1] + 2 * 0.22 * G * ds));
    for (i = M - 2; i >= 0; i--) v[i] = Math.min(v[i], Math.sqrt(v[i + 1] * v[i + 1] + 2 * 0.45 * G * ds));
    v[M - 1] = 2;
    for (i = M - 40; i < M; i++) v[i] = Math.max(1, v[i] * (M - i) / 40);
    // time along path
    var ts = [0]; for (i = 1; i < M; i++) ts.push(ts[i - 1] + ds / Math.max(0.5, (v[i] + v[i - 1]) / 2));
    var Ttot = ts[M - 1], H = 0.02, NN = Math.floor(Ttot / H);
    var sAt = [], j = 0;
    for (var q = 0; q < NN; q++) { var t = q * H; while (j < M - 2 && ts[j + 1] < t) j++; var f = (t - ts[j]) / (ts[j + 1] - ts[j]); sAt.push(j + f); }
    function at(arr, s) { var a = Math.floor(s), b = Math.min(M - 1, a + 1), f = s - a; return arr[a] + (arr[b] - arr[a]) * f; }
    var vt = sAt.map(function (s) { return at(v, s); });
    var lonT = deriv(vt, H).map(function (x) { return x / G; });
    // driver smoothness varies by segment: rough segments get less filtering
    var rough = sAt.map(function (s) { return (Math.floor(s / 420) % 3 === 1) ? 1 : 0; });
    var lonS = smooth(lonT, 25), lonR = smooth(lonT.map(function (x) { return x < 0 ? x * 1.5 : x * 1.3; }), 2);
    lonT = lonT.map(function (x, k) { return rough[k] ? lonR[k] : lonS[k]; });
    var latT = sAt.map(function (s, k) { return vt[k] * vt[k] * at(kap, s) / G; });
    var latR = latT.map(function (x, k) { return x + (rough[k] ? 0.12 * Math.sin(k * H * 5) * Math.min(1, Math.abs(x) * 3) : 0); });
    // GPS at 1 Hz
    var lat0 = opts.lat0 || 5.5, lon0 = opts.lon0 || 101.4, kx = 111320 * Math.cos(lat0 * Math.PI / 180), ky = 110540;
    var gps = [], imu = [];
    for (q = 0; q < NN; q += 50) {
      var s0 = sAt[q];
      var xx = at(px, s0) + gauss(r) * 2, yy = at(py, s0) + gauss(r) * 2;
      gps.push([Math.round(q * H * 1000 + 400), +(lat0 + yy / ky).toFixed(6), +(lon0 + xx / kx).toFixed(6), +Math.max(0, vt[q] + gauss(r) * 0.15).toFixed(2), null, 4, 300 + 0.03 * s0]);
    }
    // IMU 50 Hz with arbitrary mount rotation (phone upright, yawed)
    var yawM = 0.6, pitch = 1.2;
    for (q = 0; q < NN; q++) {
      var fwd = lonT[q] * G, lt = latR[q] * G, up = G;
      // car frame (f, l, u) -> phone frame
      var a1 = fwd * Math.cos(yawM) - lt * Math.sin(yawM), a2 = fwd * Math.sin(yawM) + lt * Math.cos(yawM);
      var b1 = a1, b2 = a2 * Math.cos(pitch) - up * Math.sin(pitch), b3 = a2 * Math.sin(pitch) + up * Math.cos(pitch);
      imu.push([Math.round(q * H * 1000), +(b1 + gauss(r) * 0.35).toFixed(3), +(b2 + gauss(r) * 0.35).toFixed(3), +(b3 + gauss(r) * 0.6).toFixed(3)]);
    }
    return { format: 'heeltoe-rec', v: 1, startedAt: opts.startedAt || Date.UTC(2026, 10, 22, 0, 30), device: 'example', gps: gps, imu: imu };
  }

  return { process: process, parseFile: parseFile, parseRecording: parseRecording, parseGPX: parseGPX, insights: insights, pack: pack, unpack: unpack, slimCorner: slimCorner, synth: synth, SUB: SUB, score: score };
})();
if (typeof module !== 'undefined') module.exports = HT;
