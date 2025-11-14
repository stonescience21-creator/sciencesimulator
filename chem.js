/* chem.js — Beast Mode chemistry engine + EVERYTHING (solubility, redox, ionic eqns, stability, multi-step, conditions, builder extras)
   - Conceptual only. No lab protocols, no temperatures.
   - Designed to plug into chem.html from previous steps (periodic builder present).
*/

(() => {
  // ---------------- CONFIG / DATA ----------------
  // Common reagents and polyatomic ions
  const COMMON_OPTIONS = [
    "H2O","HCl","NaOH","NaCl","H2SO4","H2","O2","N2","CO2","CH4","C2H6",
    "CuSO4","AgNO3","KI","KBr","BaCl2","Pb(NO3)2","Na2CO3","CaCO3",
    "Fe","Zn","Cu","Ag","Na","K","NH3","HNO3","KOH","Mg","Al","C","S"
  ];

  const POLYATOMIC = [
    {key:"SO4", formula:"SO4", charge:-2},
    {key:"NO3", formula:"NO3", charge:-1},
    {key:"NH4", formula:"NH4", charge:1},
    {key:"CO3", formula:"CO3", charge:-2},
    {key:"PO4", formula:"PO4", charge:-3},
    {key:"OH", formula:"OH", charge:-1}
  ];

  // Solubility rules (simplified but expanded)
  // We'll encode rules as functions that say whether a salt (cation, anion) is soluble.
  const SOLUBILITY_RULES = {
    // nitrates always soluble
    nitrate: (cation, anion) => anion.includes("NO3"),
    // group 1 cations always soluble
    group1: (cation, anion) => ["Li","Na","K","Rb","Cs"].includes(cation),
    // ammonium soluble
    ammonium: (cation, anion) => cation === "NH4",
    // chlorides bromides iodides generally soluble except Ag, Pb, Hg
    halides: (cation, anion) => {
      if (!["Cl","Br","I"].some(x => anion.includes(x))) return false;
      return !["Ag","Pb","Hg"].includes(cation);
    },
    // sulfates mostly soluble except Ba, Sr, Pb, Ca (simplified)
    sulfate: (cation, anion) => {
      if (!anion.includes("SO4")) return false;
      return !["Ba","Sr","Pb","Ca"].includes(cation);
    },
    // carbonates, phosphates: insoluble except group1 and NH4
    carbonate_phosphate: (cation, anion) => {
      if (!anion.includes("CO3") && !anion.includes("PO4")) return true;
      return !["Li","Na","K","Rb","Cs","NH4"].includes(cation);
    },
    // hydroxides: soluble for group1; heavier alkaline earth slightly soluble (Ca,Sr,Ba trickle)
    hydroxide: (cation, anion) => {
      if (!anion.includes("OH")) return true;
      if (["Li","Na","K","Rb","Cs","NH4"].includes(cation)) return true;
      if (["Ca","Sr","Ba"].includes(cation)) return false; // treat as insoluble for precipitation decisions
      return false;
    }
  };

  // Common oxidation numbers and valences (naive)
  const COMMON_OX = {
    H: 1, Li:1, Na:1, K:1, Rb:1, Cs:1,
    Be:2, Mg:2, Ca:2, Sr:2, Ba:2,
    O: -2, F: -1, Cl: -1, Br: -1, I: -1,
    N: -3, S: -2, P: -3,
    Al: 3, Zn: 2, Ag: 1, Fe: 2, Fe3:3, Cu: 2, Pb: 2
  };

  // Periodic grid (compact) — used by builder tooltip
  const ELEMENT_INFO = {
    H:{name:"Hydrogen",z:1,group:"nonmetal",valence:1},
    He:{name:"Helium",z:2,group:"noble",valence:0},
    Li:{name:"Lithium",z:3,group:"alkali",valence:1},
    Be:{name:"Beryllium",z:4,group:"alkaline",valence:2},
    B:{name:"Boron",z:5,group:"metalloid",valence:3},
    C:{name:"Carbon",z:6,group:"nonmetal",valence:4},
    N:{name:"Nitrogen",z:7,group:"nonmetal",valence:3},
    O:{name:"Oxygen",z:8,group:"nonmetal",valence:2},
    F:{name:"Fluorine",z:9,group:"halogen",valence:1},
    Ne:{name:"Neon",z:10,group:"noble",valence:0},
    Na:{name:"Sodium",z:11,group:"alkali",valence:1},
    Mg:{name:"Magnesium",z:12,group:"alkaline",valence:2},
    Al:{name:"Aluminium",z:13,group:"post-transition",valence:3},
    Si:{name:"Silicon",z:14,group:"metalloid",valence:4},
    P:{name:"Phosphorus",z:15,group:"nonmetal",valence:3},
    S:{name:"Sulfur",z:16,group:"nonmetal",valence:2},
    Cl:{name:"Chlorine",z:17,group:"halogen",valence:1},
    Ar:{name:"Argon",z:18,group:"noble",valence:0},
    K:{name:"Potassium",z:19,group:"alkali",valence:1},
    Ca:{name:"Calcium",z:20,group:"alkaline",valence:2},
    Fe:{name:"Iron",z:26,group:"transition",valence:2},
    Cu:{name:"Copper",z:29,group:"transition",valence:2},
    Zn:{name:"Zinc",z:30,group:"transition",valence:2},
    Ag:{name:"Silver",z:47,group:"transition",valence:1},
    Pb:{name:"Lead",z:82,group:"post-transition",valence:2}
    // expand if needed
  };

  // Reaction conditions heuristics database
  const REACTION_CONDITIONS_DB = {
    "Combustion": ["heat (ignition)", "oxygen"],
    "Acid-Base Neutralization": ["aqueous solution"],
    "Precipitation": ["aqueous solution"],
    "Single Displacement": ["metal reactivity", "aqueous solution"],
    "No notable reaction": []
  };

  // ---------------- UTILITIES ----------------

  function toSubscript(eq) {
    // numbers -> sub, simple charge patterns -> superscript
    return eq
      .replace(/(\d+)/g, "<sub>$1</sub>")
      .replace(/(\^?\-?\+?\d*[\+\-]$)|([+\-]\d+$)/g, (m)=>`<sup>${m.replace('^','')}</sup>`);
  }

  function stripCharges(formula) {
    // Remove ending charge markers like 2+, 3-, +, - and any caret '^'
    return (formula || "").replace(/\^/g,"").replace(/(\d+[\+\-]$)|([\+\-]$)|([\+\-]\d+$)/g,"").trim();
  }

  // parse formula with parentheses: returns element counts {C:2,H:6,...}
  function parseFormula(formula) {
    formula = (formula||"").trim();
    formula = stripCharges(formula);
    if (!formula) return {};
    let i = 0;
    function parseSegment() {
      const counts = {};
      while (i < formula.length) {
        const ch = formula[i];
        if (ch === "(") {
          i++;
          const inner = parseSegment();
          if (formula[i] === ")") i++;
          const mult = parseNumber() || 1;
          multiply(inner, mult);
          merge(counts, inner);
        } else if (ch === ")") {
          break;
        } else {
          const el = parseElement();
          if (!el) break;
          const num = parseNumber() || 1;
          counts[el] = (counts[el] || 0) + num;
        }
      }
      return counts;
    }
    function parseElement() {
      if (!/[A-Z]/.test(formula[i])) return null;
      let el = formula[i++];
      if (i < formula.length && /[a-z]/.test(formula[i])) el += formula[i++];
      return el;
    }
    function parseNumber() {
      let s = "";
      while (i < formula.length && /\d/.test(formula[i])) s += formula[i++];
      return s ? parseInt(s,10) : null;
    }
    function multiply(obj, m) { for (const k in obj) obj[k] = obj[k] * m; }
    function merge(a,b) { for (const k in b) a[k] = (a[k]||0) + b[k]; }
    return parseSegment();
  }

  function gcd(a,b){ return b ? gcd(b,a%b) : Math.abs(a); }

  // Balancer: numeric RREF-ish approach (same general idea as before)
  function balanceEquation(reactants, products) {
    const species = reactants.concat(products);
    const parsed = species.map(s => parseFormula(stripCharges(s)));
    const elements = Array.from(new Set(parsed.flatMap(p => Object.keys(p))));
    if (elements.length === 0) return null;
    const mat = elements.map(el => parsed.map((p,idx)=> idx<reactants.length ? (p[el]||0) : -(p[el]||0)));

    // RREF
    const A = mat.map(r => r.slice());
    const rows = A.length, cols = A[0].length;
    const EPS = 1e-9;
    let r = 0, pivotCols = [];
    for (let c=0; c<cols && r<rows; c++) {
      let pivot = r;
      while (pivot<rows && Math.abs(A[pivot][c]) < EPS) pivot++;
      if (pivot === rows) continue;
      [A[r],A[pivot]] = [A[pivot],A[r]];
      const pv = A[r][c];
      for (let j=c;j<cols;j++) A[r][j] /= pv;
      for (let i2=0;i2<rows;i2++) if (i2!==r) {
        const f = A[i2][c];
        if (Math.abs(f)>EPS) for (let j=c;j<cols;j++) A[i2][j] -= f*A[r][j];
      }
      pivotCols.push(c); r++;
    }
    const solution = new Array(cols).fill(0);
    const free = [];
    for (let c=0;c<cols;c++) if (!pivotCols.includes(c)) free.push(c);
    if (free.length === 0) return null;
    free.forEach(f=>solution[f]=1);
    for (let i=pivotCols.length-1;i>=0;i--) {
      const pc = pivotCols[i];
      const rowIdx = A.findIndex(row => Math.abs(row[pc]-1) < EPS);
      if (rowIdx === -1) continue;
      let val = 0;
      for (let j=pc+1;j<cols;j++) val += A[rowIdx][j]*solution[j];
      solution[pc] = -val;
    }
    // scale to integer
    let factor = 1;
    for (let s=1;s<=80;s++){
      const scaled = solution.map(x=>x*s);
      if (scaled.every(v=>Math.abs(v-Math.round(v))<1e-6)) { factor = s; break; }
    }
    let scaled = solution.map(x=>Math.round(x*factor));
    const sign = scaled.find(v=>v>0)?1:-1;
    scaled = scaled.map(v=>v*sign);
    const g = scaled.reduce((a,b)=>gcd(a,b));
    return scaled.map(v=>v/g);
  }

  // ---------------- SOLUBILITY ENGINE ----------------
  function isSoluble(cation, anion) {
    // check rules — if any rule returns true, treat soluble; if a rule says insoluble, treat insoluble
    // We'll run specific checks and return boolean
    // Nitrates always soluble
    if (SOLUBILITY_RULES.nitrate(null, anion)) return true;
    if (SOLUBILITY_RULES.group1(cation, anion)) return true;
    if (SOLUBILITY_RULES.ammonium(cation, anion)) return true;
    if (SOLUBILITY_RULES.halides(cation, anion)) return true;
    if (SOLUBILITY_RULES.sulfate(cation, anion)) return true;
    // carbonates/phosphates mostly insoluble
    if (!SOLUBILITY_RULES.carbonate_phosphate(cation, anion)) return false;
    // hydroxide logic: if OH present and not group1, tend to insoluble for precipitation
    if (!SOLUBILITY_RULES.hydroxide(cation, anion)) return false;

    // fallback: consider soluble
    return true;
  }

  // Given two ionic formulas (like AgNO3 and NaCl), attempt to compute products and decide precipitate
  function precipitationProducts(a,b) {
    const A = splitIonic(a);
    const B = splitIonic(b);
    if (!A || !B) return null;
    const prod1 = A.cation + B.anion;
    const prod2 = B.cation + A.anion;
    const insol1 = !isSoluble(A.cation, B.anion);
    const insol2 = !isSoluble(B.cation, A.anion);
    const products = [];
    if (insol1) products.push(prod1);
    if (insol2) products.push(prod2);
    return products.length ? products : null;
  }

  // Ionic split helper
  function splitIonic(formula) {
    // naive: split first element symbol as cation; remainder as anion (strip numbers)
    const f = stripCharges(formula);
    const m = f.match(/^([A-Z][a-z]?)(.*)$/);
    if (!m) return null;
    const cat = m[1];
    const an = (m[2]||"").replace(/\d+/g,"").replace(/^\(|\)$/g,"");
    return {cation:cat, anion:an||""};
  }

  // ---------------- REDOX MODULE ----------------
  // Heuristic oxidation state estimation (naive)
  function estimateOxStates(atomCounts) {
    const elems = Object.keys(atomCounts);
    const ox = {};
    let sumKnown = 0;
    let unknowns = [];

    elems.forEach(e => {
      if (e === "O") { ox[e] = -2; sumKnown += ox[e]*atomCounts[e]; }
      else if (e === "H") { ox[e] = 1; sumKnown += ox[e]*atomCounts[e]; }
      else if (COMMON_OX[e] !== undefined) { ox[e] = COMMON_OX[e]; sumKnown += ox[e]*atomCounts[e]; }
      else unknowns.push(e);
    });

    if (unknowns.length === 1) {
      const u = unknowns[0];
      ox[u] = Math.round((-sumKnown/atomCounts[u])*10)/10;
    } else unknowns.forEach(u=>ox[u]=0);

    return ox;
  }

  function detectRedoxChange(reactants, products) {
    // reactants/products: arrays of formula strings
    // For each species, estimate oxidation states and compare elements that change
    const parsedReact = reactants.map(r => ({f:r, p:parseFormula(stripCharges(r)), ox:estimateOxStates(parseFormula(stripCharges(r)))}));
    const parsedProd  = products.map(r => ({f:r, p:parseFormula(stripCharges(r)), ox:estimateOxStates(parseFormula(stripCharges(r)))}));

    // sum oxidation numbers of each element across reactants and products
    const sumOxReact = {};
    parsedReact.forEach(spec => {
      for (const e in spec.p) sumOxReact[e] = (sumOxReact[e]||0) + (spec.ox[e]||0)*spec.p[e];
    });
    const sumOxProd = {};
    parsedProd.forEach(spec => {
      for (const e in spec.p) sumOxProd[e] = (sumOxProd[e]||0) + (spec.ox[e]||0)*spec.p[e];
    });

    const changes = [];
    for (const e of new Set([...Object.keys(sumOxReact), ...Object.keys(sumOxProd)])) {
      const rVal = sumOxReact[e] || 0;
      const pVal = sumOxProd[e] || 0;
      if (Math.abs(rVal - pVal) > 1e-6) {
        changes.push({element:e, reactSum:rVal, prodSum:pVal, delta:pVal - rVal});
      }
    }

    // classify which elements lost vs gained electrons: delta <0 means decreased oxidation sum (reduction of total?), careful
    // We'll interpret per-atom change by comparing average oxidation per atom if counts differ — approximation
    return {changes, parsedReact, parsedProd};
  }

  // Build half-reactions heuristic (very naive): find element(s) with biggest change and propose left/right electron transfer
  function buildHalfReactions(reactants, products) {
    const det = detectRedoxChange(reactants, products);
    if (!det.changes.length) return null;
    // pick top two elements by absolute delta
    const sorted = det.changes.sort((a,b)=>Math.abs(b.delta)-Math.abs(a.delta));
    // create textual half reactions
    const halves = sorted.slice(0,2).map(ch => {
      const sign = ch.delta>0 ? "oxidized (lost electrons)" : "reduced (gained electrons)";
      return `${ch.element}: net change ${ch.delta.toFixed(1)} → ${sign}`;
    });
    return {summary: halves.join("; "), details: det};
  }

  // ---------------- IONIC & NET IONIC EQUATIONS ----------------
  function ionicRepresentation(formula) {
    // naive split into ions for simple ionic salts like NaCl, AgNO3, CuSO4
    // returns array of ions as strings, e.g., ["Na+","Cl-"] or null if cannot split
    // Try to find cation and anion in simple format
    const s = stripCharges(formula);
    const m = s.match(/^([A-Z][a-z]?)(.*)$/);
    if (!m) return null;
    const c = m[1];
    let an = (m[2]||"").replace(/\d+/g,"");
    if (!an) return null;
    // assign charges crudely: cation charge = common ox for element or +1; anion charge deduced to balance (coarse)
    const cOx = COMMON_OX[c] || 1;
    // find anion charge by summing known oxidation states in polyatomic (very naive)
    const parsedAn = parseFormula(an);
    let anCharge = 0;
    for (const el in parsedAn) {
      const ox = (COMMON_OX[el] !== undefined) ? COMMON_OX[el] : 0;
      anCharge += ox * parsedAn[el];
    }
    // total compound neutral: cOx + anCharge = 0 => anChargeNeg = -cOx
    // this simplistic method often fails, but provides a rough ionic labeling
    const cion = c + (cOx>1?cOx:"+");
    const anionLabel = an + (anCharge? (anCharge>0?("+"+anCharge):anCharge) : "-");
    return [cion, anionLabel];
  }

  function spectatorCancellation(ionsLeft, ionsRight) {
    // ionsLeft and ionsRight arrays of ions like "Na+", "Cl-"
    // Count occurrences and cancel common spectators to produce net ions
    const count = (arr)=> {
      const map = {};
      arr.forEach(i=> map[i] = (map[i]||0)+1);
      return map;
    };
    const L = count(ionsLeft), R = count(ionsRight);
    const netLeft = [], netRight = [];
    // cancel matching keys
    const keys = new Set([...Object.keys(L), ...Object.keys(R)]);
    keys.forEach(k=>{
      const l = L[k]||0; const r = R[k]||0;
      if (l>r) {
        for (let i=0;i<l-r;i++) netLeft.push(k);
      } else if (r>l) {
        for (let i=0;i<r-l;i++) netRight.push(k);
      }
    });
    return {netLeft, netRight};
  }

  // ---------------- STABILITY CHECKER ----------------
  function isPlausibleCompound(formula) {
    // checks for impossible stoichiometry w.r.t common valences (very heuristic)
    const parsed = parseFormula(stripCharges(formula));
    if (!Object.keys(parsed).length) return false;
    // If only single element and count>0, it's fine (elemental)
    if (Object.keys(parsed).length === 1) return true;
    // attempt to compute a rough charge balance using common oxidation states
    let total = 0, atoms = 0;
    for (const el in parsed) {
      const count = parsed[el];
      const ox = COMMON_OX[el] !== undefined ? COMMON_OX[el] : 0;
      total += ox * count;
      atoms += count;
    }
    // stable ionic compound should have total near 0; for covalent, accept nonzero
    if (Math.abs(total) <= 4) return true;
    // flag improbable if total huge
    return false;
  }

  // ---------------- HEAT / ENERGY ESTIMATOR ----------------
  function estimateHeatOfReaction(type, reactants, products) {
    // conceptual mapping only
    const map = {
      "Combustion": {rating:5, note:"strongly exothermic"},
      "Acid-Base Neutralization": {rating:4, note:"exothermic"},
      "Precipitation": {rating:2, note:"slightly exothermic"},
      "Single Displacement": {rating:3, note:"moderately exothermic or endothermic depending on metals"},
      "No notable reaction": {rating:1, note:"no significant heat change predicted"}
    };
    return map[type] || {rating:1, note:"no prediction"};
  }

  // ---------------- MULTI-STEP SIMULATOR ----------------
  function multiStepSimulate(initialReagents, maxSteps=4) {
    // Applies predictor iteratively: products become new reagents if they are reactive with something
    let steps = [];
    let current = initialReagents.slice();
    for (let s=0;s<maxSteps;s++) {
      const pred = predictReaction(current);
      steps.push({step:s+1, reagents:current.slice(), prediction:pred});
      if (!pred || !pred.equation || pred.type === "No notable reaction") break;
      // parse products from balanced equation if present
      const prods = extractProductsFromEquation(pred.equation);
      // if prods are same as reagents or empty, break
      if (!prods.length) break;
      // next set: use products that are not identical to any current reagent
      current = prods;
    }
    return steps;
  }

  function extractProductsFromEquation(equation) {
    // naive parsing of "A + B → C + D"
    if (!equation) return [];
    const parts = equation.split("→");
    if (parts.length < 2) return [];
    const right = parts[1].trim();
    const prods = right.split("+").map(x=>x.trim()).map(x=>x.replace(/^\d+\s*/,"").trim());
    return prods;
  }

  // ---------------- REACTION PREDICTOR (integrated upgrades) ----------------
  function finalBalance(react, prod) {
    const c = balanceEquation(react, prod);
    if (!c) return null;
    return formatEquation(react, prod, c);
  }

  function formatEquation(react, prod, coeffs) {
    const left = react.map((r,i) => ((coeffs[i] && coeffs[i] !== 1) ? coeffs[i] + " " : "") + r).join(" + ");
    const right = prod.map((p,i) => ((coeffs[react.length + i] && coeffs[react.length + i] !== 1) ? coeffs[react.length + i] + " " : "") + p).join(" + ");
    return left + " → " + right;
  }

  function predictReaction(reagents) {
    // Input: array of strings possibly with charges; we'll use both raw and stripped.
    const raw = reagents.slice();
    const clean = reagents.map(r => stripCharges(r));
    // Combustion
    const hasO2 = clean.includes("O2");
    const hydrocarbon = clean.find(r => {
      const p = parseFormula(r);
      return (p.C || 0) > 0 && (p.H || 0) > 0;
    });
    if (hydrocarbon && hasO2) {
      const react = [hydrocarbon, "O2"];
      const prod = ["CO2","H2O"];
      const eq = finalBalance(react, prod);
      const heat = estimateHeatOfReaction("Combustion", react, prod);
      return {
        type:"Combustion",
        equation:eq,
        explanation:`Combustion of ${hydrocarbon} in oxygen produces CO2 and H2O (conceptual).`,
        heat, conditions:REACTION_CONDITIONS_DB["Combustion"]
      };
    }

    // Acid-Base (more robust): detect H+ donor and OH- donor
    const acids = clean.filter(r => /^H[A-Z0-9(]/.test(r) || ["HCl","H2SO4","HNO3"].includes(r));
    const bases = clean.filter(r => /OH/.test(r) || ["NaOH","KOH","Ca(OH)2"].includes(r));
    if (acids.length && bases.length) {
      const acid = acids[0], base = bases[0];
      // try to form salt: naive anion = acid without H; cation = base cation
      const anion = acid.replace(/^H/,"") || "X";
      const cation = (base.match(/^([A-Z][a-z]?)/) || ["M"])[0];
      const salt = cation + anion;
      const react = [acid, base];
      const prod = [salt, "H2O"];
      const eq = finalBalance(react, prod);
      const ionicLeft = ionicRepresentation(acid) || [acid];
      const ionicRight = ionicRepresentation(base) || [base];
      const redox = buildHalfReactions(react, prod);
      const heat = estimateHeatOfReaction("Acid-Base Neutralization", react, prod);
      return {
        type:"Acid-Base Neutralization",
        equation:eq,
        explanation:`Naive neutralization: ${salt} + H2O (conceptual).`,
        ionic:{left:ionicLeft,right:ionicRight},
        netIonic:computeNetIonic(react, prod),
        redox,
        heat,
        conditions:REACTION_CONDITIONS_DB["Acid-Base Neutralization"]
      };
    }

    // Precipitation: check all pairs for insoluble combos (improved solubility rules)
    for (let i=0;i<clean.length;i++) {
      for (let j=i+1;j<clean.length;j++) {
        const prods = precipitationProducts(clean[i], clean[j]);
        if (prods) {
          const eq = finalBalance([clean[i],clean[j]], prods);
          const net = computeNetIonic([clean[i],clean[j]], prods);
          const heat = estimateHeatOfReaction("Precipitation",[clean[i],clean[j]],prods);
          return {
            type:"Precipitation",
            equation:eq,
            explanation:`Possible precipitation forming ${prods.join(", ")} using solubility rules.`,
            products:prods,
            netIonic: net,
            heat,
            conditions:REACTION_CONDITIONS_DB["Precipitation"]
          };
        }
      }
    }

    // Single displacement (metal reactivity list)
    const metals = ["K","Na","Ca","Mg","Al","Zn","Fe","Pb","Cu","Ag"];
    for (let i=0;i<clean.length;i++) {
      for (let j=0;j<clean.length;j++) if (i!==j) {
        const a = clean[i], b = clean[j];
        if (metals.includes(a) && splitIonic(b) && metals.includes(splitIonic(b).cation)) {
          const salt = splitIonic(b);
          if (metals.indexOf(a) < metals.indexOf(salt.cation)) {
            const newSalt = a + salt.anion;
            const eq = finalBalance([a,b],[newSalt, salt.cation]);
            const heat = estimateHeatOfReaction("Single Displacement",[a,b],[newSalt, salt.cation]);
            return { type:"Single Displacement", equation:eq, explanation:`${a} may displace ${salt.cation} to form ${newSalt}.`, heat, conditions:REACTION_CONDITIONS_DB["Single Displacement"] };
          }
        }
      }
    }

    // Synthesis / Decomposition heuristics: detect if two elements combine to form one compound
    if (clean.length === 2) {
      // try to form a single product by concatenation heuristic if plausible
      const prod = [clean.join("")];
      if (isPlausibleCompound(prod[0])) {
        const eq = finalBalance(clean, prod);
        if (eq) return { type:"Synthesis (heuristic)", equation:eq, explanation:"Simple combination predicted (naive).", conditions:["heat maybe"] };
      }
    }
    if (clean.length === 1) {
      // attempt decomposition into simple products like metal oxide -> metal + O2 (very naive)
      const single = clean[0];
      const p = parseFormula(single);
      if (p && p.O && Object.keys(p).length>1) {
        // propose decomposition: compound -> oxide fragments (not robust)
        // skip for safety if we can't produce reasonable products
      }
    }

    // Nothing flagged
    return { type:"No notable reaction", equation:null, explanation:"No reaction predicted by the expanded engine.", conditions:[] };
  }

  // ---------------- NET IONIC COMPUTATION ----------------
  function computeNetIonic(reactants, products) {
    // try to represent reactants and products as ions (very naive)
    // Return {completeIonicLeft, completeIonicRight, netIonicLeft, netIonicRight}
    const leftIons = reactants.flatMap(r => ionicSpeciesFor(r) || [r]);
    const rightIons = products.flatMap(r => ionicSpeciesFor(r) || [r]);
    const cancellation = spectatorCancellation(leftIons, rightIons);
    return {completeLeft:leftIons, completeRight:rightIons, netLeft:cancellation.netLeft, netRight:cancellation.netRight};
  }

  function ionicSpeciesFor(formula) {
    // For salts of the form MAn (very naive)
    const sp = splitIonic(formula);
    if (!sp) return null;
    // choose charges heuristically
    const catOx = COMMON_OX[sp.cation] || 1;
    // approximate anion charge by summing common ox of its elements
    const parsedAn = parseFormula(sp.anion);
    let anCharge = 0;
    for (const el in parsedAn) {
      const count = parsedAn[el];
      const ox = COMMON_OX[el] !== undefined ? COMMON_OX[el] : 0;
      anCharge += ox*count;
    }
    // compute integer charges to balance (coarse)
    let cCharge = catOx;
    // return strings with signs
    const cLabel = sp.cation + (cCharge===1?"+": (cCharge>1?String(cCharge)+"+":"+"));
    const aLabel = sp.anion + (anCharge? (anCharge>0?("+"+anCharge): (""+anCharge)) : "-");
    // return array; for a polyatomic with counts you'd normally also show stoichiometry — skip for simplicity
    return [cLabel, aLabel];
  }

  // ---------------- BUILDER EXTRAS & UI ----------------
  // DOM elements (assumes chem.html layout with these IDs)
  const reagentArea = document.getElementById("reagent-area");
  const addBtn = document.getElementById("add-reagent");
  const removeBtn = document.getElementById("remove-reagent");
  const predictBtn = document.getElementById("predict");
  const resultDiv = document.getElementById("result");

  // Builder popup elements (periodic builder)
  const popup = document.getElementById("periodic-popup");
  const periodicTableContainer = document.getElementById("periodic-table");
  const polyListDiv = document.getElementById("poly-list");
  const rawDiv = document.getElementById("formula-raw");
  const prettyDiv = document.getElementById("formula-pretty");
  const insertBtn = document.getElementById("insert-formula");
  const cancelBtn = document.getElementById("cancel-formula");
  const closeBtn = document.getElementById("periodic-close");

  // Builder internal state
  let builder = { parts: [], cursorSlot: null };
  function resetBuilder() { builder.parts = []; updateBuilderDisplay(); }

  function addElementToBuilder(sym) {
    const last = builder.parts[builder.parts.length-1];
    if (last && last.type==="el" && last.sym === sym) {
      last.count = (last.count||1) + 1;
    } else {
      builder.parts.push({type:"el", sym, count:1});
    }
    updateBuilderDisplay();
  }
  function addPolyatomic(formula) {
    builder.parts.push({type:"poly", formula, count:1});
    updateBuilderDisplay();
  }
  function addParenGroup() {
    builder.parts.push({type:"group", parts:[], count:1});
    updateBuilderDisplay();
  }
  function addSubscriptToLast() {
    const last = builder.parts[builder.parts.length-1]; if (!last) return;
    last.count = (last.count||1) + 1; updateBuilderDisplay();
  }
  function addChargeToLast() {
    const last = builder.parts[builder.parts.length-1]; if (!last) return;
    last.charge = last.charge ? (last.charge==="+"?"2+":(last.charge==="2?"?"-":"+")) : "+";
    updateBuilderDisplay();
  }

  function renderRaw() {
    // produce linear raw formula string
    return builder.parts.map(p=>{
      if (p.type==="el") return p.sym + (p.count && p.count>1 ? p.count : "") + (p.charge?p.charge:"");
      if (p.type==="poly") return p.formula + (p.count && p.count>1 ? p.count : "");
      if (p.type==="group") {
        const inner = (p.parts && p.parts.length) ? p.parts.map(ip => ip.sym ? ip.sym + (ip.count>1?ip.count:"") : ip.formula).join("") : "";
        return "(" + inner + ")" + (p.count && p.count>1 ? p.count : "");
      }
      return "";
    }).join("");
  }
  function updateBuilderDisplay() {
    const raw = renderRaw();
    rawDiv.textContent = raw || "(empty)";
    prettyDiv.innerHTML = toSubscript(raw || "");
  }

  // Build periodic grid
  const PERIODIC = [
    [{sym:"H"},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{sym:"He"}],
    [{sym:"Li"},{sym:"Be"},{},{},{},{},{},{},{sym:"B"},{sym:"C"},{sym:"N"},{sym:"O"},{sym:"F"},{sym:"Ne"},{},{},{},{}],
    [{sym:"Na"},{sym:"Mg"},{},{},{},{},{},{},{sym:"Al"},{sym:"Si"},{sym:"P"},{sym:"S"},{sym:"Cl"},{sym:"Ar"},{},{},{},{}],
    [{sym:"K"},{sym:"Ca"},{sym:"Sc"},{sym:"Ti"},{sym:"V"},{sym:"Cr"},{sym:"Mn"},{sym:"Fe"},{sym:"Co"},{sym:"Ni"},{sym:"Cu"},{sym:"Zn"},{sym:"Ga"},{sym:"Ge"},{sym:"As"},{sym:"Se"},{sym:"Br"},{sym:"Kr"}]
  ];

  function buildPeriodicGrid() {
    if (!periodicTableContainer) return;
    periodicTableContainer.innerHTML = "";
    PERIODIC.forEach(row => {
      row.forEach(cell => {
        const tile = document.createElement("div");
        tile.className = "periodic-cell small";
        if (!cell || !cell.sym) { tile.style.visibility = "hidden"; periodicTableContainer.appendChild(tile); return; }
        tile.textContent = cell.sym;
        tile.title = `${cell.sym} — ${ELEMENT_INFO[cell.sym]?.name || ""}`;
        tile.onclick = () => addElementToBuilder(cell.sym);
        tile.onmouseover = () => showElementTooltip(cell.sym);
        tile.onmouseleave = hideElementTooltip;
        periodicTableContainer.appendChild(tile);
      });
    });
  }

  // poly list
  function buildPolyList() {
    if (!polyListDiv) return;
    polyListDiv.innerHTML = "";
    POLYATOMIC.forEach(p=>{
      const b = document.createElement("button");
      b.className = "btn-ghost small";
      b.textContent = p.key;
      b.onclick = () => addPolyatomic(p.formula);
      polyListDiv.appendChild(b);
    });
  }

  // element tooltip small overlay (basic)
  let tooltipEl = null;
  function showElementTooltip(sym) {
    if (!ELEMENT_INFO[sym]) return;
    hideElementTooltip();
    tooltipEl = document.createElement("div");
    tooltipEl.style.position = "fixed";
    tooltipEl.style.background = "#131416";
    tooltipEl.style.border = "1px solid rgba(255,255,255,0.04)";
    tooltipEl.style.padding = "6px 8px";
    tooltipEl.style.borderRadius = "6px";
    tooltipEl.style.color = "#ddd";
    tooltipEl.style.zIndex = 2000;
    tooltipEl.textContent = `${sym} — ${ELEMENT_INFO[sym].name} (Z=${ELEMENT_INFO[sym].z}) valence ${ELEMENT_INFO[sym].valence}`;
    document.body.appendChild(tooltipEl);
    // place near mouse? too late — center top
    tooltipEl.style.left = (window.event ? (window.event.clientX + 10) : 200) + "px";
    tooltipEl.style.top = (window.event ? (window.event.clientY + 10) : 100) + "px";
  }
  function hideElementTooltip() { if (tooltipEl) { tooltipEl.remove(); tooltipEl = null; } }

  function initPopupControls() {
    document.querySelectorAll(".popup-inner [data-action]").forEach(btn=>{
      btn.addEventListener("click", (e)=>{
        const act = e.currentTarget.dataset.action;
        if (act === "clear") resetBuilder();
        else if (act === "open-paren") addParenGroup();
        else if (act === "add-sub") addSubscriptToLast();
        else if (act === "charge-plus") addChargeToLast();
        updateBuilderDisplay();
      });
    });
    cancelBtn.onclick = () => { popup.classList.add("hidden"); resetBuilder(); };
    closeBtn.onclick = cancelBtn.onclick;
  }

  function openPeriodicBuilder(slotIndex, onInsert) {
    builder.cursorSlot = slotIndex;
    popup.classList.remove("hidden");
    resetBuilder();
    buildPeriodicGrid();
    buildPolyList();
    initPopupControls();
    insertBtn.onclick = () => {
      const formula = rawDiv.textContent === "(empty)" ? "" : rawDiv.textContent;
      popup.classList.add("hidden");
      if (onInsert) onInsert(formula);
      resetBuilder();
    };
  }

  // ---------------- UI SLOTS ----------------
  let slotCount = 4;
  const MAX_SLOTS = 6;
  const MIN_SLOTS = 1;

  function buildSlots() {
    if (!reagentArea) return;
    reagentArea.innerHTML = "";
    for (let s=0;s<slotCount;s++) {
      const div = document.createElement("div");
      div.className = "reagent";
      div.innerHTML = `
        <label>Reagent ${s+1}</label>
        <div style="display:flex;gap:6px;align-items:center;">
          <input class="formula-input" placeholder="Formula or use builder..." value="${COMMON_OPTIONS[s]||''}" style="flex:1">
          <button class="btn-ghost small open-builder" data-slot="${s}">Builder</button>
        </div>
        <div class="mini">Parsed: <span class="parsed">-</span></div>
      `;
      reagentArea.appendChild(div);

      const input = div.querySelector(".formula-input");
      const parsed = div.querySelector(".parsed");
      input.addEventListener("input", ()=> {
        try {
          const p = parseFormula(input.value||"");
          parsed.textContent = Object.entries(p).map(([k,v])=>k+(v>1?v:"")).join(", ") || "-";
        } catch { parsed.textContent = "-"; }
      });
      const b = div.querySelector(".open-builder");
      b.onclick = () => openPeriodicBuilder(s, (formula) => {
        if (!formula) return;
        input.value = formula;
        input.dispatchEvent(new Event('input'));
      });
      input.dispatchEvent(new Event('input'));
    }
  }

  addBtn.onclick = () => {
    if (slotCount >= MAX_SLOTS) return alert("Max reagents: 6");
    slotCount++; buildSlots();
  };
  removeBtn.onclick = () => {
    if (slotCount <= MIN_SLOTS) return alert("At least one reagent required.");
    slotCount--; buildSlots();
  };

  // ---------------- PREDICT BUTTON + OUTPUT FORMATTING ----------------
  predictBtn.onclick = () => {
    const inputs = Array.from(document.querySelectorAll(".formula-input")).map(i=>i.value.trim()).filter(v=>v);
    if (!inputs.length) { resultDiv.innerHTML = "<em>No reagents added.</em>"; return; }
    try {
      // Multi-step simulation included: present single-step prediction plus multi-step summary
      const singlePred = predictReaction(inputs);
      const multi = multiStepSimulate(inputs, 5);
      // build output HTML
      let out = "";
      out += `<div style="margin-bottom:12px;"><strong>Type:</strong> ${singlePred.type}</div>`;
      if (singlePred.equation) {
        out += `<div style="margin-bottom:12px;"><strong>Equation:</strong><br>${singlePred.equation}<br><br><strong>Formatted:</strong><br>${toSubscript(singlePred.equation)}</div>`;
      }
      out += `<div style="margin-bottom:12px;"><strong>Explanation:</strong><br>${singlePred.explanation}</div>`;
      // conditions
      if (singlePred.conditions && singlePred.conditions.length) {
        out += `<div style="margin-bottom:12px;"><strong>Conditions / Catalysts:</strong><br>${singlePred.conditions.join(", ")}</div>`;
      }
      // heat
      if (singlePred.heat) {
        out += `<div style="margin-bottom:12px;"><strong>Conceptual Heat:</strong> ${"★".repeat(singlePred.heat.rating)} (${singlePred.heat.note})</div>`;
      }
      // ionic / net ionic
      if (singlePred.netIonic) {
        out += `<div style="margin-bottom:12px;"><strong>Net Ionic (approx):</strong><br>`;
        const ni = singlePred.netIonic;
        out += `<em>Complete Ionic Left:</em> ${ni.completeLeft.join(", ")}<br>`;
        out += `<em>Complete Ionic Right:</em> ${ni.completeRight.join(", ")}<br>`;
        out += `<em>Net Left:</em> ${ni.netLeft.join(", ")}<br><em>Net Right:</em> ${ni.netRight.join(", ")}<br></div>`;
      } else if (singlePred.ionic) {
        out += `<div style="margin-bottom:12px;"><strong>Approx Ionic:</strong><br>Left: ${singlePred.ionic.left.join(", ")}<br>Right: ${singlePred.ionic.right.join(", ")}</div>`;
      }

      // redox summary
      if (singlePred.redox) {
        out += `<div style="margin-bottom:12px;"><strong>Redox summary (heuristic):</strong><br>${singlePred.redox.summary || ""}</div>`;
      }

      // parsed reagents and oxidation states
      out += `<div style="margin-bottom:8px;"><strong>Parsed reagents & estimated oxidation states:</strong></div>`;
      inputs.forEach(f => {
        const parsed = parseFormula(stripCharges(f));
        const ox = estimateOxStates(parsed);
        out += `<div style="margin-bottom:6px;"><strong>${escapeHtml(f)}</strong><br>Elements: ${Object.entries(parsed).map(([k,v])=>k+(v>1?v:"")).join(", ") || "(none)"}<br>Oxidation (est): ${Object.entries(ox).map(([k,v])=>k+":"+v).join(", ")}</div>`;
      });

      // multi-step summary if >1 step
      if (multi && multi.length>1) {
        out += `<div style="margin-top:12px;"><strong>Multi-step simulation:</strong><br>`;
        multi.forEach(step=>{
          out += `<div style="margin-top:6px;"><em>Step ${step.step}:</em> Reagents: ${step.reagents.join(" + ")}<br>`;
          out += `Predicted: ${step.prediction.type}${step.prediction.equation?(" — "+step.prediction.equation):""}</div>`;
        });
        out += `</div>`;
      }

      // stability checks
      out += `<div style="margin-top:12px;"><strong>Stability checks:</strong><br>`;
      inputs.forEach(f=>{
        const ok = isPlausibleCompound(f);
        out += `${escapeHtml(f)}: ${ok?'<span style="color:#9fe59f">Plausible</span>':'<span style="color:#f39a9a">Possibly implausible</span>'}<br>`;
      });
      out += `</div>`;

      resultDiv.innerHTML = out;
    } catch (e) {
      resultDiv.innerHTML = "<strong>Error:</strong> " + escapeHtml(e.message || String(e));
    }
  };

  // escape helper
  function escapeHtml(s) {
    return String(s).replace(/[&<>"']/g, (m)=>({ '&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;' })[m]);
  }

  // ---------------- INIT ----------------
  buildPeriodicGrid();
  buildPolyList();
  buildSlots();

  // small exposure for debugging in console (if you want)
  window.__chem = {
    parseFormula, balanceEquation, predictReaction, multiStepSimulate, isPlausibleCompound
  };
resultDiv.classList.remove("reaction-flash");
void resultDiv.offsetWidth; // force reflow
resultDiv.classList.add("reaction-flash");

// ---------------- TUTORIAL ----------------
(function() {
  const tutorialSteps = [
    {
      text: "Welcome to the Chemistry Simulator! Let's learn the basics.",
      highlight: null
    },
    {
      text: "These boxes are where you enter reagents.",
      highlight: () => document.querySelector(".reagent")
    },
    {
      text: "Use the 'Builder' button to construct formulas using the periodic table.",
      highlight: () => document.querySelector(".open-builder")
    },
    {
      text: "Add more reagents with this button.",
      highlight: () => document.getElementById("add-reagent")
    },
    {
      text: "Remove excess reagent slots here.",
      highlight: () => document.getElementById("remove-reagent")
    },
    {
      text: "Press Predict to simulate your reaction.",
      highlight: () => document.getElementById("predict")
    },
    {
      text: "And you're done! Go experiment.",
      highlight: null
    }
  ];

  let tutIndex = 0;
  const tutOverlay = document.getElementById("tut-overlay");
  const tutText = document.getElementById("tut-text");
  const tutNext = document.getElementById("tut-next");

  function runTutorial() {
    tutOverlay.classList.remove("hidden");
    showTutorialStep();
  }

  function showTutorialStep() {
    document.querySelectorAll(".tut-highlight").forEach(el => el.classList.remove("tut-highlight"));
    const step = tutorialSteps[tutIndex];
    tutText.innerHTML = step.text;

    if (step.highlight) {
      const target = step.highlight();
      if (target) target.classList.add("tut-highlight");
    }
  }

  tutNext.onclick = () => {
    tutIndex++;
    if (tutIndex >= tutorialSteps.length) {
      tutOverlay.classList.add("hidden");
      localStorage.setItem("chemTutorialDone","1");
      document.querySelectorAll(".tut-highlight").forEach(el => el.classList.remove("tut-highlight"));
      return;
    }
    showTutorialStep();
  };

  // Run once unless already completed
  if (!localStorage.getItem("chemTutorialDone")) {
    setTimeout(runTutorial, 700);
  }

  // EXPOSE TO GLOBAL SCOPE
  window.openChemTutorial = () => {
    tutIndex = 0;
    runTutorial();
  };
})();

})();
