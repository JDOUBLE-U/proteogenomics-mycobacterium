/*
 *  /MathJax/jax/output/HTML-CSS/autoload/mmultiscripts.js
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
    a.mmultiscripts.Augment({
        toHTML: function (L, J, E) {
            L = this.HTMLcreateSpan(L);
            var R = this.HTMLgetScale();
            var n = b.createStack(L), g;
            var l = b.createBox(n);
            if (this.data[this.base]) {
                var m = this.data[this.base].toHTML(l);
                if (E != null) {
                    this.data[this.base].HTMLstretchV(l, J, E)
                } else {
                    if (J != null) {
                        this.data[this.base].HTMLstretchH(l, J)
                    }
                }
                b.Measured(m, l)
            } else {
                l.bbox = this.HTMLzeroBBox()
            }
            var P = b.TeX.x_height * R, C = b.TeX.scriptspace * R * 0.75;
            var B = this.HTMLgetScripts(n, C);
            var o = B[0], f = B[1], y = B[2], k = B[3];
            var j = R;
            for (var M = 1; M < this.data.length; M++) {
                if (this.data[M] && this.data[M].spanID) {
                    j = this.data[M].HTMLgetScale();
                    break
                }
            }
            var G = b.TeX.sup_drop * j, F = b.TeX.sub_drop * j;
            var z = l.bbox.h - G, x = l.bbox.d + F, Q = 0, I;
            if (l.bbox.ic) {
                Q = l.bbox.ic
            }
            if (this.data[this.base] && (this.data[this.base].type === "mi" || this.data[this.base].type === "mo")) {
                if (this.data[this.base].data.join("").length === 1 && l.bbox.scale === 1 && !this.data[this.base].Get("largeop")) {
                    z = x = 0
                }
            }
            var K = this.getValues("subscriptshift", "superscriptshift"), H = this.HTMLgetMu(L);
            K.subscriptshift = (K.subscriptshift === "" ? 0 : b.length2em(K.subscriptshift, H));
            K.superscriptshift = (K.superscriptshift === "" ? 0 : b.length2em(K.superscriptshift, H));
            var w = 0;
            if (y) {
                w = y.bbox.w + Q
            } else {
                if (k) {
                    w = k.bbox.w - Q
                }
            }
            if (w < 0) {
                w = 0
            }
            b.placeBox(l, w, 0);
            if (!f && !k) {
                x = Math.max(x, b.TeX.sub1 * R, K.subscriptshift);
                if (o) {
                    x = Math.max(x, o.bbox.h - (4 / 5) * P)
                }
                if (y) {
                    x = Math.max(x, y.bbox.h - (4 / 5) * P)
                }
                if (o) {
                    b.placeBox(o, w + l.bbox.w + C - Q, -x)
                }
                if (y) {
                    b.placeBox(y, 0, -x)
                }
            } else {
                if (!o && !y) {
                    g = this.getValues("displaystyle", "texprimestyle");
                    I = b.TeX[(g.displaystyle ? "sup1" : (g.texprimestyle ? "sup3" : "sup2"))];
                    z = Math.max(z, I * R, K.superscriptshift);
                    if (f) {
                        z = Math.max(z, f.bbox.d + (1 / 4) * P)
                    }
                    if (k) {
                        z = Math.max(z, k.bbox.d + (1 / 4) * P)
                    }
                    if (f) {
                        b.placeBox(f, w + l.bbox.w + C, z)
                    }
                    if (k) {
                        b.placeBox(k, 0, z)
                    }
                } else {
                    x = Math.max(x, b.TeX.sub2 * R);
                    var A = b.TeX.rule_thickness * R;
                    var N = (o || y).bbox.h, O = (f || k).bbox.d;
                    if (y) {
                        N = Math.max(N, y.bbox.h)
                    }
                    if (k) {
                        O = Math.max(O, k.bbox.d)
                    }
                    if ((z - O) - (N - x) < 3 * A) {
                        x = 3 * A - z + O + N;
                        G = (4 / 5) * P - (z - O);
                        if (G > 0) {
                            z += G;
                            x -= G
                        }
                    }
                    z = Math.max(z, K.superscriptshift);
                    x = Math.max(x, K.subscriptshift);
                    if (f) {
                        b.placeBox(f, w + l.bbox.w + C, z)
                    }
                    if (k) {
                        b.placeBox(k, w + Q - k.bbox.w, z)
                    }
                    if (o) {
                        b.placeBox(o, w + l.bbox.w + C - Q, -x)
                    }
                    if (y) {
                        b.placeBox(y, w - y.bbox.w, -x)
                    }
                }
            }
            this.HTMLhandleSpace(L);
            this.HTMLhandleColor(L);
            var e = L.bbox;
            e.dx = w;
            e.s = C;
            e.u = z;
            e.v = x;
            e.delta = Q;
            e.px = w + l.bbox.w;
            return L
        }, HTMLgetScripts: function (p, q) {
            var o, d, e = [];
            var n = 1, g = this.data.length, f = 0;
            for (var h = 0; h < 4; h += 2) {
                while (n < g && this.data[n].type !== "mprescripts") {
                    for (var l = h; l < h + 2; l++) {
                        if (this.data[n] && this.data[n].type !== "none") {
                            if (!e[l]) {
                                e[l] = b.createBox(p);
                                e[l].bbox = this.HTMLemptyBBox({});
                                if (f) {
                                    b.createBlank(e[l], f);
                                    e[l].bbox.w = e[l].bbox.rw = f
                                }
                            }
                            this.data[n].toHTML(e[l]);
                            this.HTMLcombineBBoxes(this.data[n], e[l].bbox)
                        }
                        n++
                    }
                    d = e[h];
                    o = e[h + 1];
                    if (d && o) {
                        if (d.bbox.w < o.bbox.w) {
                            b.createBlank(d, o.bbox.w - d.bbox.w);
                            d.bbox.w = o.bbox.w;
                            d.bbox.rw = Math.max(d.bbox.w, d.bbox.rw)
                        } else {
                            if (d.bbox.w > o.bbox.w) {
                                b.createBlank(o, d.bbox.w - o.bbox.w);
                                o.bbox.w = d.bbox.w;
                                o.bbox.rw = Math.max(o.bbox.w, o.bbox.rw)
                            }
                        }
                    }
                    if (d) {
                        f = d.bbox.w
                    } else {
                        if (o) {
                            f = o.bbox.w
                        }
                    }
                }
                n++;
                f = 0
            }
            for (l = 0; l < 4; l++) {
                if (e[l]) {
                    e[l].bbox.w += q;
                    e[l].bbox.rw = Math.max(e[l].bbox.w, e[l].bbox.rw);
                    e[l].bbox.name = (["sub", "sup", "presub", "presup"])[l];
                    this.HTMLcleanBBox(e[l].bbox)
                }
            }
            return e
        }, HTMLstretchH: a.mbase.HTMLstretchH, HTMLstretchV: a.mbase.HTMLstretchV
    });
    MathJax.Hub.Startup.signal.Post("HTML-CSS mmultiscripts Ready");
    MathJax.Ajax.loadComplete(b.autoloadDir + "/mmultiscripts.js")
});