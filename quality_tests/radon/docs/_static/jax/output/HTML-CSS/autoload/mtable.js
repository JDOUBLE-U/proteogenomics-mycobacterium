/*
 *  /MathJax/jax/output/HTML-CSS/autoload/mtable.js
 *
 *  Copyright (c) 2009-2014 The MathJax Consortium
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 */

MathJax.Hub.Register.StartupHook("HTML-CSS Jax Ready", function () {
    var c = "2.4.0";
    var a = MathJax.ElementJax.mml, b = MathJax.OutputJax["HTML-CSS"];
    a.mtable.Augment({
        toHTML: function (t) {
            t = this.HTMLcreateSpan(t);
            if (this.data.length === 0) {
                return t
            }
            var K = this.getValues("columnalign", "rowalign", "columnspacing", "rowspacing", "columnwidth", "equalcolumns", "equalrows", "columnlines", "rowlines", "frame", "framespacing", "align", "useHeight", "width", "side", "minlabelspacing");
            var aI = K.width.match(/%$/);
            var ax = b.createStack(t);
            var aF = this.HTMLgetScale(), az = this.HTMLgetMu(t), aA = -1;
            var ap = [], at = [], ak = [], av = [], au = [], af, ae, ao = -1, ad, an, Z, aE, R, aB, aN = [], aS;
            var I = b.FONTDATA.lineH * aF * K.useHeight, O = b.FONTDATA.lineD * aF * K.useHeight;
            for (af = 0, ad = this.data.length; af < ad; af++) {
                aE = this.data[af];
                Z = (aE.type === "mlabeledtr" ? aA : 0);
                av[af] = [];
                ap[af] = I;
                at[af] = O;
                for (ae = Z, an = aE.data.length + Z; ae < an; ae++) {
                    if (ak[ae] == null) {
                        if (ae > ao) {
                            ao = ae
                        }
                        au[ae] = b.createStack(b.createBox(ax));
                        ak[ae] = -b.BIGDIMEN
                    }
                    av[af][ae] = b.createBox(au[ae]);
                    aN.push(aE.data[ae - Z].toHTML(av[af][ae]))
                }
            }
            b.MeasureSpans(aN);
            for (af = 0, ad = this.data.length; af < ad; af++) {
                aE = this.data[af];
                Z = (aE.type === "mlabeledtr" ? aA : 0);
                for (ae = Z, an = aE.data.length + Z; ae < an; ae++) {
                    R = aE.data[ae - Z];
                    if (R.isMultiline) {
                        av[af][ae].style.width = "100%"
                    }
                    if (R.isEmbellished()) {
                        aB = R.CoreMO();
                        var aR = aB.Get("minsize", true);
                        if (aR) {
                            var aK = aB.HTMLspanElement().bbox;
                            if (aB.HTMLcanStretch("Vertical")) {
                                aS = aK.h + aK.d;
                                if (aS) {
                                    aR = b.length2em(aR, az, aS);
                                    if (aR * aK.h / aS > ap[af]) {
                                        ap[af] = aR * aK.h / aS
                                    }
                                    if (aR * aK.d / aS > at[af]) {
                                        at[af] = aR * aK.d / aS
                                    }
                                }
                            } else {
                                if (aB.HTMLcanStretch("Horizontal")) {
                                    aR = b.length2em(aR, az, aK.w);
                                    if (aR > ak[ae]) {
                                        ak[ae] = aR
                                    }
                                }
                            }
                        }
                    }
                    if (av[af][ae].bbox.h > ap[af]) {
                        ap[af] = av[af][ae].bbox.h
                    }
                    if (av[af][ae].bbox.d > at[af]) {
                        at[af] = av[af][ae].bbox.d
                    }
                    if (av[af][ae].bbox.w > ak[ae]) {
                        ak[ae] = av[af][ae].bbox.w
                    }
                }
            }
            var aD = MathJax.Hub.SplitList;
            var ay = aD(K.columnspacing), aP = aD(K.rowspacing), e = aD(K.columnalign), E = aD(K.rowalign), d = aD(K.columnlines), z = aD(K.rowlines), aL = aD(K.columnwidth), V = [];
            for (af = 0, ad = ay.length; af < ad; af++) {
                ay[af] = b.length2em(ay[af], az)
            }
            for (af = 0, ad = aP.length; af < ad; af++) {
                aP[af] = b.length2em(aP[af], az)
            }
            while (ay.length < ao) {
                ay.push(ay[ay.length - 1])
            }
            while (e.length <= ao) {
                e.push(e[e.length - 1])
            }
            while (d.length < ao) {
                d.push(d[d.length - 1])
            }
            while (aL.length <= ao) {
                aL.push(aL[aL.length - 1])
            }
            while (aP.length < av.length) {
                aP.push(aP[aP.length - 1])
            }
            while (E.length <= av.length) {
                E.push(E[E.length - 1])
            }
            while (z.length < av.length) {
                z.push(z[z.length - 1])
            }
            if (au[aA]) {
                e[aA] = (K.side.substr(0, 1) === "l" ? "left" : "right");
                ay[aA] = -ak[aA]
            }
            for (af = 0, ad = av.length; af < ad; af++) {
                aE = this.data[af];
                V[af] = [];
                if (aE.rowalign) {
                    E[af] = aE.rowalign
                }
                if (aE.columnalign) {
                    V[af] = aD(aE.columnalign);
                    while (V[af].length <= ao) {
                        V[af].push(V[af][V[af].length - 1])
                    }
                }
            }
            if (K.equalrows) {
                var aC = Math.max.apply(Math, ap), X = Math.max.apply(Math, at);
                for (af = 0, ad = av.length; af < ad; af++) {
                    Z = ((aC + X) - (ap[af] + at[af])) / 2;
                    ap[af] += Z;
                    at[af] += Z
                }
            }
            aS = ap[0] + at[av.length - 1];
            for (af = 0, ad = av.length - 1; af < ad; af++) {
                aS += Math.max(0, at[af] + ap[af + 1] + aP[af])
            }
            var aH = 0, aG = 0, aU, g = aS;
            if (K.frame !== "none" || (K.columnlines + K.rowlines).match(/solid|dashed/)) {
                var w = aD(K.framespacing);
                if (w.length != 2) {
                    w = aD(this.defaults.framespacing)
                }
                aH = b.length2em(w[0], az);
                aG = b.length2em(w[1], az);
                g = aS + 2 * aG
            }
            var aj, aT, ab = "";
            if (typeof(K.align) !== "string") {
                K.align = String(K.align)
            }
            if (K.align.match(/(top|bottom|center|baseline|axis)( +(-?\d+))?/)) {
                ab = RegExp.$3;
                K.align = RegExp.$1
            } else {
                K.align = this.defaults.align
            }
            if (ab !== "") {
                ab = parseInt(ab);
                if (ab < 0) {
                    ab = av.length + 1 + ab
                }
                if (ab < 1) {
                    ab = 1
                } else {
                    if (ab > av.length) {
                        ab = av.length
                    }
                }
                aj = 0;
                aT = -(aS + aG) + ap[0];
                for (af = 0, ad = ab - 1; af < ad; af++) {
                    var N = Math.max(0, at[af] + ap[af + 1] + aP[af]);
                    aj += N;
                    aT += N
                }
            } else {
                aj = ({
                    top: -(ap[0] + aG),
                    bottom: aS + aG - ap[0],
                    center: aS / 2 - ap[0],
                    baseline: aS / 2 - ap[0],
                    axis: aS / 2 + b.TeX.axis_height * aF - ap[0]
                })[K.align];
                aT = ({
                    top: -(aS + 2 * aG),
                    bottom: 0,
                    center: -(aS / 2 + aG),
                    baseline: -(aS / 2 + aG),
                    axis: b.TeX.axis_height * aF - aS / 2 - aG
                })[K.align]
            }
            var ac, ag = 0, B = 0, L = 0, aa = 0, ah = 0, am = [], ar = [], S = 1;
            if (K.equalcolumns && K.width !== "auto") {
                if (aI) {
                    ac = (100 / (ao + 1)).toFixed(2).replace(/\.?0+$/, "") + "%";
                    for (af = 0, ad = Math.min(ao + 1, aL.length); af < ad; af++) {
                        aL[af] = ac
                    }
                    ac = 0;
                    ag = 1;
                    ah = ao + 1;
                    for (af = 0, ad = Math.min(ao + 1, ay.length); af < ad; af++) {
                        ac += ay[af]
                    }
                } else {
                    ac = b.length2em(K.width, az);
                    for (af = 0, ad = Math.min(ao + 1, ay.length); af < ad; af++) {
                        ac -= ay[af]
                    }
                    ac /= ao + 1;
                    for (af = 0, ad = Math.min(ao + 1, aL.length); af < ad; af++) {
                        ak[af] = ac
                    }
                }
            } else {
                for (af = 0, ad = Math.min(ao + 1, aL.length); af < ad; af++) {
                    if (aL[af] === "auto") {
                        B += ak[af]
                    } else {
                        if (aL[af] === "fit") {
                            ar[ah] = af;
                            ah++;
                            B += ak[af]
                        } else {
                            if (aL[af].match(/%$/)) {
                                am[aa] = af;
                                aa++;
                                L += ak[af];
                                ag += b.length2em(aL[af], az, 1)
                            } else {
                                ak[af] = b.length2em(aL[af], az);
                                B += ak[af]
                            }
                        }
                    }
                }
                if (aI) {
                    ac = 0;
                    for (af = 0, ad = Math.min(ao, ay.length); af < ad; af++) {
                        ac += ay[af]
                    }
                    if (ag > 0.98) {
                        S = 0.98 / ag;
                        ag = 0.98
                    }
                } else {
                    if (K.width === "auto") {
                        if (ag > 0.98) {
                            S = L / (B + L);
                            ac = B + L
                        } else {
                            ac = B / (1 - ag)
                        }
                    } else {
                        ac = b.length2em(K.width, az);
                        for (af = 0, ad = Math.min(ao + 1, ay.length); af < ad; af++) {
                            ac -= ay[af]
                        }
                    }
                    for (af = 0, ad = am.length; af < ad; af++) {
                        ak[am[af]] = b.length2em(aL[am[af]], az, ac * S);
                        B += ak[am[af]]
                    }
                    if (Math.abs(ac - B) > 0.01) {
                        if (ah && ac > B) {
                            ac = (ac - B) / ah;
                            for (af = 0, ad = ar.length; af < ad; af++) {
                                ak[ar[af]] += ac
                            }
                        } else {
                            ac = ac / B;
                            for (ae = 0; ae <= ao; ae++) {
                                ak[ae] *= ac
                            }
                        }
                    }
                    if (K.equalcolumns) {
                        var Q = Math.max.apply(Math, ak);
                        for (ae = 0; ae <= ao; ae++) {
                            ak[ae] = Q
                        }
                    }
                }
            }
            var T = aj, o, r, aQ;
            Z = (au[aA] ? aA : 0);
            for (ae = Z; ae <= ao; ae++) {
                for (af = 0, ad = av.length; af < ad; af++) {
                    if (av[af][ae]) {
                        Z = (this.data[af].type === "mlabeledtr" ? aA : 0);
                        R = this.data[af].data[ae - Z];
                        if (R.HTMLcanStretch("Horizontal")) {
                            av[af][ae].bbox = R.HTMLstretchH(au[ae], ak[ae]).bbox
                        } else {
                            if (R.HTMLcanStretch("Vertical")) {
                                aB = R.CoreMO();
                                var aJ = aB.symmetric;
                                aB.symmetric = false;
                                av[af][ae].bbox = R.HTMLstretchV(au[ae], ap[af], at[af]).bbox;
                                av[af][ae].HH = null;
                                if (av[af][ae].bbox.h > ap[af]) {
                                    av[af][ae].bbox.H = av[af][ae].bbox.h;
                                    av[af][ae].bbox.h = ap[af]
                                }
                                if (av[af][ae].bbox.d > at[af]) {
                                    av[af][ae].bbox.D = av[af][ae].bbox.d;
                                    av[af][ae].bbox.d = at[af]
                                }
                                aB.symmetric = aJ
                            }
                        }
                        aQ = R.rowalign || this.data[af].rowalign || E[af];
                        o = ({
                            top: ap[af] - av[af][ae].bbox.h,
                            bottom: av[af][ae].bbox.d - at[af],
                            center: ((ap[af] - at[af]) - (av[af][ae].bbox.h - av[af][ae].bbox.d)) / 2,
                            baseline: 0,
                            axis: 0
                        })[aQ] || 0;
                        aQ = (R.columnalign || V[af][ae] || e[ae]);
                        b.alignBox(av[af][ae], aQ, T + o)
                    }
                    if (af < av.length - 1) {
                        T -= Math.max(0, at[af] + ap[af + 1] + aP[af])
                    }
                }
                T = aj
            }
            if (aI) {
                var G = b.createBox(ax);
                G.style.left = G.style.top = 0;
                G.style.right = b.Em(ac + 2 * aH);
                G.style.display = "inline-block";
                G.style.height = "0px";
                if (b.msieRelativeWidthBug) {
                    G = b.createBox(G);
                    G.style.position = "relative";
                    G.style.height = "1em";
                    G.style.width = "100%";
                    G.bbox = ax.bbox
                }
                var aO = 0, aV = aH, k, l;
                if (ah) {
                    k = 100 * (1 - ag) / ah, l = B / ah
                } else {
                    k = 100 * (1 - ag) / (ao + 1);
                    l = B / (ao + 1)
                }
                for (ae = 0; ae <= ao; ae++) {
                    b.placeBox(au[ae].parentNode, 0, 0);
                    au[ae].style.position = "relative";
                    au[ae].style.left = b.Em(aV);
                    au[ae].style.width = "100%";
                    au[ae].parentNode.parentNode.removeChild(au[ae].parentNode);
                    var al = b.createBox(G);
                    b.addBox(al, au[ae]);
                    au[ae] = al;
                    var h = al.style;
                    h.display = "inline-block";
                    h.left = aO + "%";
                    if (aL[ae].match(/%$/)) {
                        var u = parseFloat(aL[ae]) * S;
                        if (ah === 0) {
                            h.width = (k + u) + "%";
                            aO += k + u;
                            al = b.createBox(al);
                            b.addBox(al, au[ae].firstChild);
                            al.style.left = 0;
                            al.style.right = b.Em(l);
                            aV -= l
                        } else {
                            h.width = u + "%";
                            aO += u
                        }
                    } else {
                        if (aL[ae] === "fit" || ah === 0) {
                            h.width = k + "%";
                            al = b.createBox(al);
                            b.addBox(al, au[ae].firstChild);
                            al.style.left = 0;
                            al.style.right = b.Em(l - ak[ae]);
                            aV += ak[ae] - l;
                            aO += k
                        } else {
                            h.width = b.Em(ak[ae]);
                            aV += ak[ae]
                        }
                    }
                    if (b.msieRelativeWidthBug) {
                        b.addText(al.firstChild, b.NBSP);
                        al.firstChild.style.position = "relative"
                    }
                    aV += ay[ae];
                    if (d[ae] !== "none" && ae < ao && ae !== aA) {
                        r = b.createBox(G);
                        r.style.left = aO + "%";
                        r = b.createRule(r, g, 0, 1.25 / b.em);
                        r.style.position = "absolute";
                        r.bbox = {h: g, d: 0, w: 0, rw: 1.25 / b.em, lw: 0};
                        r.parentNode.bbox = ax.bbox;
                        b.placeBox(r, aV - ay[ae] / 2, aT, true);
                        r.style.borderStyle = d[ae]
                    }
                }
            } else {
                var U = aH;
                for (ae = 0; ae <= ao; ae++) {
                    if (!au[ae].bbox.width) {
                        b.setStackWidth(au[ae], ak[ae])
                    }
                    if (aL[ae] !== "auto" && aL[ae] !== "fit") {
                        au[ae].bbox.width = ak[ae];
                        au[ae].bbox.isFixed = true
                    }
                    b.placeBox(au[ae].parentNode, U, 0);
                    U += ak[ae] + ay[ae];
                    if (d[ae] !== "none" && ae < ao && ae !== aA) {
                        r = b.createRule(ax, g, 0, 1.25 / b.em);
                        b.addBox(ax, r);
                        r.bbox = {h: g, d: 0, w: 0, rw: 1.25 / b.em, lw: 0};
                        b.placeBox(r, U - ay[ae] / 2, aT, true);
                        r.style.borderStyle = d[ae]
                    }
                }
            }
            ax.bbox.d = -aT;
            ax.bbox.h = g + aT;
            b.setStackWidth(ax, ax.bbox.w + aH);
            aU = ax.bbox.w;
            var ai;
            if (K.frame !== "none") {
                ai = b.createFrame(ax, g, 0, aU, 1.25 / b.em, K.frame);
                b.addBox(ax, ai);
                b.placeBox(ai, 0, aT, true);
                if (aI) {
                    ai.style.width = "100%"
                }
            }
            T = aj;
            for (af = 0, ad = av.length - 1; af < ad; af++) {
                o = Math.max(0, at[af] + ap[af + 1] + aP[af]);
                if (z[af] !== "none") {
                    r = b.createRule(ax, 1.25 / b.em, 0, aU);
                    b.addBox(ax, r);
                    r.bbox = {h: 1.25 / b.em, d: 0, w: aU, rw: aU, lw: 0};
                    b.placeBox(r, 0, T - at[af] - (o - at[af] - ap[af + 1]) / 2, true);
                    if (z[af] === "dashed" || aI) {
                        r.style.borderTop = r.style.height + " " + z[af];
                        r.style.height = 0;
                        r.style.width = r.style.borderLeftWidth;
                        r.style.borderLeft = "";
                        if (aI) {
                            r.style.width = "100%"
                        }
                    }
                }
                T -= o
            }
            if (aI) {
                t.bbox.width = K.width;
                ax.style.width = "100%"
            }
            if (au[aA]) {
                var aw = ax.bbox.w, q;
                var aq = this.getValues("indentalignfirst", "indentshiftfirst", "indentalign", "indentshift");
                if (aq.indentalignfirst !== a.INDENTALIGN.INDENTALIGN) {
                    aq.indentalign = aq.indentalignfirst
                }
                if (aq.indentalign === a.INDENTALIGN.AUTO) {
                    aq.indentalign = this.displayAlign
                }
                if (aq.indentshiftfirst !== a.INDENTSHIFT.INDENTSHIFT) {
                    aq.indentshift = aq.indentshiftfirst
                }
                if (aq.indentshift === "auto") {
                    aq.indentshift = this.displayIndent
                }
                var aM = b.createStack(t, false, "100%");
                b.addBox(aM, ax);
                b.alignBox(ax, aq.indentalign, 0);
                if (aq.indentshift && aq.indentalign !== a.INDENTALIGN.CENTER) {
                    q = b.length2em(aq.indentshift, az);
                    aw += q;
                    ax.style[aq.indentalign] = b.Em(q)
                }
                au[aA].parentNode.parentNode.removeChild(au[aA].parentNode);
                b.addBox(aM, au[aA]);
                b.alignBox(au[aA], e[aA], 0);
                if (b.msieRelativeWidthBug) {
                    ax.style.top = au[aA].style.top = ""
                }
                if (aI) {
                    ax.style.width = K.width;
                    t.bbox.width = "100%"
                }
                q = b.length2em(K.minlabelspacing, az);
                au[aA].style.marginRight = au[aA].style.marginLeft = b.Em(q);
                if (aq.indentalign === a.INDENTALIGN.CENTER) {
                    aw += 4 * q + 2 * au[aA].bbox.w
                } else {
                    if (aq.indentalign !== e[aA]) {
                        aw += 2 * q + au[aA].bbox.w
                    }
                }
                t.style.minWidth = t.bbox.minWidth = aM.style.minWidth = aM.bbox.minWidth = b.Em(aw)
            }
            if (!aI) {
                this.HTMLhandleSpace(t)
            }
            var v = this.HTMLhandleColor(t);
            if (v && aI) {
                if (!ai) {
                    ai = b.createFrame(ax, g, 0, aU, 0, "none");
                    b.addBox(ax, ai);
                    b.placeBox(ai, 0, aT, true);
                    ai.style.width = "100%"
                }
                ai.style.backgroundColor = v.style.backgroundColor;
                ai.parentNode.insertBefore(ai, ai.parentNode.firstChild);
                v.parentNode.removeChild(v)
            }
            return t
        }, HTMLhandleSpace: function (d) {
            d.bbox.keepPadding = true;
            d.bbox.exact = true;
            if (!this.hasFrame && d.bbox.width == null) {
                d.style.paddingLeft = d.style.paddingRight = b.Em(1 / 6)
            }
            this.SUPER(arguments).HTMLhandleSpace.call(this, d)
        }
    });
    a.mtd.Augment({
        toHTML: function (e, d, g) {
            e = this.HTMLcreateSpan(e);
            if (this.data[0]) {
                var f = this.data[0].toHTML(e);
                if (g != null) {
                    f = this.data[0].HTMLstretchV(e, d, g)
                } else {
                    if (d != null) {
                        f = this.data[0].HTMLstretchH(e, d)
                    }
                }
                e.bbox = f.bbox
            }
            this.HTMLhandleSpace(e);
            this.HTMLhandleColor(e);
            return e
        }, HTMLstretchH: a.mbase.HTMLstretchH, HTMLstretchV: a.mbase.HTMLstretchV
    });
    MathJax.Hub.Startup.signal.Post("HTML-CSS mtable Ready");
    MathJax.Ajax.loadComplete(b.autoloadDir + "/mtable.js")
});