"""Run with python -m cocycle; no native extension is required."""
import argparse,json
from pathlib import Path

def main():
    parser=argparse.ArgumentParser(description='Exact verification and complete mathematical case reports.')
    sub=parser.add_subparsers(dest='command',required=True)
    p=sub.add_parser('verify');p.add_argument('--output',default='results');p.add_argument('--quick',action='store_true',help='360 positive commutations and named loop examples only')
    p=sub.add_parser('replay');p.add_argument('--certificates',default='results/local-certificates.json')
    p=sub.add_parser('show');p.add_argument('family',choices=['tetra','cube','commute']);p.add_argument('--case',default='1',help='tetra: s; cube: template,rotation; commute: first,second,arrangement');p.add_argument('--beta',default='beta1',choices=['beta1','beta2','beta3']);p.add_argument('--outside',type=int);p.add_argument('--positive',action='store_true');p.add_argument('--certificates',default='results/local-certificates.json');p.add_argument('--output',default='case.html');p.add_argument('--lang',choices=['en','zh'],default='en')
    p=sub.add_parser('loop');p.add_argument('kind',choices=['rotation','rolling','half-rolling','push','bracket','half-bracket']);p.add_argument('knot');p.add_argument('--other');p.add_argument('--output',default='loop.html');p.add_argument('--lang',choices=['en','zh'],default='en')
    p=sub.add_parser('serve');p.add_argument('--port',type=int,default=8765);p.add_argument('--results',default='results');p.add_argument('--lang',choices=['en','zh'],default='en')
    args=parser.parse_args()
    if args.command=='verify':
        from .verification import run
        from .loops import verify
        local=run(args.output,full=not args.quick);loops=verify(args.output,extended=not args.quick)
        print(json.dumps(dict(local=local,loop_cases=loops['cases'],loop_failures=loops['failures']),indent=2))
    elif args.command=='replay':
        from .verification import replay
        replay(json.loads(Path(args.certificates).read_text()));print('PASS: all exact linking-number certificates replayed.')
    elif args.command=='show':
        from .reports import local_case,commutation_case
        parts=[int(x) for x in args.case.split(',')]
        if args.family=='commute':
            if len(parts)!=3:parser.error('commute --case requires u,v,a')
            content=commutation_case(*parts,positive=args.positive,name=args.beta,lang=args.lang)
        else:
            if len(parts)!=(1 if args.family=='tetra' else 2):parser.error('wrong case format')
            bundle=json.loads(Path(args.certificates).read_text());content=local_case(bundle,args.family,parts[0] if len(parts)==1 else parts,args.beta,args.outside,lang=args.lang)
        Path(args.output).parent.mkdir(parents=True,exist_ok=True);Path(args.output).write_text(content);print(args.output)
    elif args.command=='loop':
        from .loops import calculate
        from .reports import loop_case
        r=calculate(args.kind,args.knot,args.other);Path(args.output).parent.mkdir(parents=True,exist_ok=True)
        Path(args.output).write_text(loop_case(r,lang=args.lang));Path(args.output).with_suffix('.json').write_text(json.dumps(r,indent=2));print(args.output)
        if not r['pass_formula']:raise SystemExit('FAIL: computed value disagrees with formula')
    else:
        from .server import serve
        serve(args.port,args.results,lang=args.lang)

if __name__=='__main__':main()
