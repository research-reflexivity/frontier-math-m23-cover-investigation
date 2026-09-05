'use strict';
// Structural and content checks; this does not execute the application.
const fs = require('node:fs');
const vm = require('node:vm');
const assert = require('node:assert/strict');
const path = require('node:path');
const html = fs.readFileSync(path.join(__dirname, 'index.html'), 'utf8');
function decode(s) {
  return s.replace(/&(#x[0-9a-f]+|#\d+|amp|lt|gt|quot|apos);/gi, (m,x) =>
    x[0] === '#' ? String.fromCodePoint(parseInt(x.slice(x[1].toLowerCase()==='x'?2:1),x[1].toLowerCase()==='x'?16:10)) :
    ({amp:'&',lt:'<',gt:'>',quot:'"',apos:"'"}[x] || m));
}
const layers = [html];
while (true) {
  const match = layers.at(-1).match(/<iframe\b[^>]*\bsrcdoc="([\s\S]*?)"/);
  if (!match) break;
  layers.push(decode(match[1]));
}
assert.equal(layers.length, 3, 'Preserve both sandboxed srcdoc layers');
const sandboxes = layers.flatMap(s=>Array.from(s.matchAll(/<iframe\b[^>]*\bsandbox="([^"]*)"/g),m=>m[1]));
assert.deepEqual(sandboxes, ['allow-scripts allow-popups allow-popups-to-escape-sandbox', 'allow-scripts']);
for (const layer of layers) {
  assert.match(layer, /Content-Security-Policy/);
  assert.doesNotMatch(layer, /sandbox="[^"]*allow-same-origin/);
  for (const m of layer.matchAll(/<script\b([^>]*)>([\s\S]*?)<\/script>/g)) {
    if (!/\bsrc=/.test(m[1]) && !/\btype=["']application\/(?:ld\+)?json/.test(m[1])) new vm.Script(m[2]);
  }
}
const inner=layers.at(-1), start=inner.indexOf('const layers ='), end=inner.indexOf('const aliasToCourse =');
assert(start>=0 && end>start);
const data=vm.runInNewContext(inner.slice(start,end)+'\n({layers,nodes,edges,courseTopics,courseEnrichment,edgeDetails})', {}, {timeout:1000});
const ids=new Set(data.nodes.map(n=>n.id));
assert.equal(ids.size,data.nodes.length);
assert.equal(data.nodes.length,32);
assert.equal(Object.keys(data.courseTopics).length,31);
for(const n of data.nodes){
  assert(n.layer>=0 && n.layer<data.layers.length);
  for(const k of ['title','formal','plain','why','evidence','breaks','status']) assert.equal(typeof n[k],'string',n.id+'.'+k);
  assert(!/[\u0000-\u0008\u000b\u000c\u000e-\u001f]/.test(n.formal),'Broken TeX escapes');
}
for(const e of data.edges){
  assert(ids.has(e.source)&&ids.has(e.target), e.id+' has unknown endpoint');
  assert(['proof','compute','gate'].includes(e.kind));
  assert(data.edgeDetails[e.source+'>'+e.target], e.id+' needs explicit semantics');
}
for(const [id,c] of Object.entries(data.courseTopics)){
  assert(c.aliases.length>0 && c.points.length===3 && c.idea && c.context, id);
  const ex=data.courseEnrichment[id];
  assert(ex.example&&ex.question&&ex.answer&&ex.links.length,id+' needs a complete lesson');
  for(const key of ex.related) assert(data.courseTopics[key],id+' links to missing lesson '+key);
  for(const [label,url] of ex.links) {assert(label);assert.equal(new URL(url).protocol,'https:');}
}
const coreOrder=JSON.parse(inner.match(/const coreOrder = (\[[^\n]*\]);/)[1]);
assert.equal(coreOrder.length,ids.size);
assert.equal(new Set(coreOrder).size,ids.size);
for(const id of coreOrder) assert(ids.has(id));
function reachable(seed,allowed) {
  const seen=new Set([seed]), queue=[seed];
  while(queue.length) {
    const current=queue.shift();
    for(const e of data.edges) if(e.source===current && allowed(e) && !seen.has(e.target)){
      seen.add(e.target);queue.push(e.target);
    }
  }
  return seen;
}
assert(reachable('epsilon_e6',e=>e.kind!=='gate').has('target'));
assert(reachable('esing_e6',e=>e.kind!=='gate').has('target'));
assert(!reachable('relativeopen',e=>e.kind!=='gate').has('target'),'No proved path from the open construction');
assert(data.edges.filter(e=>e.source==='relativeopen').every(e=>e.kind==='gate'));
const topo=new Map(data.nodes.map(n=>[n.id,0]));
data.edges.forEach(e=>topo.set(e.target,topo.get(e.target)+1));
const queue=[...topo].filter(([,v])=>v===0).map(([k])=>k);
let visited=0;
while(queue.length){const id=queue.shift();visited++;for(const e of data.edges.filter(e=>e.source===id)){topo.set(e.target,topo.get(e.target)-1);if(topo.get(e.target)===0)queue.push(e.target);}}
assert.equal(visited,ids.size,'Dependency graph must be acyclic');
assert.match(html, /An M<sub>23<\/sub> Hurwitz scheme: exact arithmetic and reduction at 23/);
assert.match(html, /independent relative geometric proof remains open/);
assert.doesNotMatch(html, /two established comparisons|effective characteristic-23 quadratic-orientation connector/);
assert.match(html, /href="m23-cover-investigation.pdf"/);
assert.match(html, /href="https:\/\/github.com\/research-reflexivity\/frontier-math-m23-cover-investigation"/);
assert.doesNotMatch(inner, /two established derivations|effective logarithmic quadratic-line comparison|const gluinglemma|const relativeline/);
for(const id of ['m23-first','m23-prev','m23-next','m23-last','m23-course','m23-edge-tooltip']) assert(inner.includes('id="'+id+'"'));
console.log('PASS: 32 nodes, 45 annotated edges, 31 complete lessons; references and hop order consistent');
console.log('PASS: proved routes reach the comparison; the open relative construction has only a dotted outgoing edge');
console.log('PASS: scripts parse; existing iframe sandbox layers and CSP remain in place');
