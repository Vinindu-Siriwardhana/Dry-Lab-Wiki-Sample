"""Shared pieces for build_page.py: the model map SVG and small helpers."""

MODEL_MAP = '''
<svg viewBox="0 0 920 492" class="modelmap" role="img"
     aria-label="Model map: loading feeds Model 1 release and Model 2 tissue transport,
     which split into the V14 arm (Models 3, 4, 5) and the FGF2-G3 arm (Model 6).
     The arms meet at a single junction. The V2 feedback edge is deleted.">
  <defs>
    <marker id="ar" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="6"
            markerHeight="6" orient="auto-start-reverse">
      <path d="M0,0 L10,5 L0,10 z" fill="var(--rule-strong)"/>
    </marker>
    <marker id="arj" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="6"
            markerHeight="6" orient="auto-start-reverse">
      <path d="M0,0 L10,5 L0,10 z" fill="var(--gold)"/>
    </marker>
    <marker id="arg" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="6"
            markerHeight="6" orient="auto-start-reverse">
      <path d="M0,0 L10,5 L0,10 z" fill="var(--rule)"/>
    </marker>
  </defs>

  <!-- deleted V2 feedback edge -->
  <path d="M150 405 C 40 405, 40 86, 150 86" fill="none" stroke="var(--rule)"
        stroke-width="2" stroke-dasharray="7 5" marker-end="url(#arg)"/>
  <text x="34" y="250" class="mm-del" transform="rotate(-90 34 250)">
    V2 feedback edge — DELETED
  </text>

  <!-- spine -->
  <line x1="460" y1="62" x2="460" y2="86" stroke="var(--rule-strong)" stroke-width="2" marker-end="url(#ar)"/>
  <line x1="460" y1="132" x2="460" y2="156" stroke="var(--rule-strong)" stroke-width="2" marker-end="url(#ar)"/>
  <path d="M460 202 L460 222 L250 222 L250 246" fill="none" stroke="var(--rule-strong)" stroke-width="2" marker-end="url(#ar)"/>
  <path d="M460 202 L460 222 L690 222 L690 246" fill="none" stroke="var(--rule-strong)" stroke-width="2" marker-end="url(#ar)"/>
  <line x1="250" y1="292" x2="250" y2="316" stroke="var(--rule-strong)" stroke-width="2" marker-end="url(#ar)"/>
  <line x1="250" y1="362" x2="250" y2="386" stroke="var(--rule-strong)" stroke-width="2" marker-end="url(#ar)"/>
  <line x1="690" y1="292" x2="690" y2="386" stroke="var(--rule-strong)" stroke-width="2" marker-end="url(#ar)"/>

  <!-- the junction -->
  <path d="M370 339 L600 339" fill="none" stroke="var(--gold)" stroke-width="3" marker-end="url(#arj)"/>
  <rect x="432" y="322" width="112" height="22" rx="4" fill="var(--gold-wash)"/>
  <text x="488" y="337" class="mm-junc">THE JUNCTION</text>

  <a href="#act2"><g>
    <rect x="300" y="26" width="320" height="36" rx="6" class="mm-box mm-load"/>
    <text x="460" y="49" class="mm-t">Loading &#160; C&#8339;&#8320; · C&#8355;&#8320; · L&#8341;&#8337;&#8343; · R&#8355;</text>
  </g></a>

  <a href="#act2"><g>
    <rect x="300" y="86" width="320" height="46" rx="6" class="mm-box mm-now"/>
    <text x="460" y="106" class="mm-t">Model 1 &#183; release from the gel</text>
    <text x="460" y="123" class="mm-s">predictive now</text>
  </g></a>

  <a href="#act2"><g>
    <rect x="300" y="156" width="320" height="46" rx="6" class="mm-box mm-now"/>
    <text x="460" y="176" class="mm-t">Model 2 &#183; tissue transport</text>
    <text x="460" y="193" class="mm-s">predictive now</text>
  </g></a>

  <text x="250" y="238" class="mm-arm mm-arm-v">V14 arm</text>
  <text x="690" y="238" class="mm-arm mm-arm-f">FGF2-G3 arm</text>

  <a href="#act4"><g>
    <rect x="120" y="246" width="260" height="46" rx="6" class="mm-box mm-scaled"/>
    <text x="250" y="266" class="mm-t">Model 3 &#183; MD2/TLR4 competition</text>
    <text x="250" y="283" class="mm-s">predictive in scaled form</text>
  </g></a>

  <a href="#act4"><g>
    <rect x="120" y="316" width="260" height="46" rx="6" class="mm-box mm-scaled"/>
    <text x="250" y="336" class="mm-t">Model 4 &#183; NF-&#954;B signalling</text>
    <text x="250" y="353" class="mm-s">predictive in scaled form</text>
  </g></a>

  <a href="#act4"><g>
    <rect x="120" y="386" width="260" height="50" rx="6" class="mm-box mm-data"/>
    <text x="250" y="406" class="mm-t">Model 5 &#183; iNOS &#8594; nitrite</text>
    <text x="250" y="423" class="mm-s">DATA IN HAND</text>
  </g></a>

  <a href="#act4"><g>
    <rect x="560" y="246" width="260" height="46" rx="6" class="mm-box mm-hyp"/>
    <text x="690" y="266" class="mm-t">fibroblast proliferation</text>
    <text x="690" y="283" class="mm-s">structured hypothesis</text>
  </g></a>

  <a href="#act3"><g>
    <rect x="560" y="386" width="260" height="50" rx="6" class="mm-box mm-hyp"/>
    <text x="690" y="406" class="mm-t">Model 6 &#183; Fisher&#8211;KPP closure</text>
    <text x="690" y="423" class="mm-s">UNANCHORED &#183; no planned experiment</text>
  </g></a>

  <a href="#act3"><g>
    <rect x="300" y="452" width="320" height="34" rx="6" class="mm-box mm-out"/>
    <text x="460" y="474" class="mm-t">Model 7 &#183; inverts the chain for design</text>
  </g></a>
  <path d="M250 436 L250 462 L300 462" fill="none" stroke="var(--rule-strong)" stroke-width="2"/>
  <path d="M690 436 L690 462 L620 462" fill="none" stroke="var(--rule-strong)" stroke-width="2"/>
</svg>
'''

LEGEND = '''
<ul class="maplegend">
  <li><span class="sw sw-now"></span>predictive now — every parameter sequence-derived or geometric</li>
  <li><span class="sw sw-scaled"></span>predictive in scaled form — runs before K<sub>D,P</sub> is measured</li>
  <li><span class="sw sw-data"></span>data in hand — a measured curve exists to fit</li>
  <li><span class="sw sw-hyp"></span>structured hypothesis — Tier-4 closures, nothing planned constrains them</li>
</ul>
'''
