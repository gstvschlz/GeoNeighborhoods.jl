#
# Renders every examples/*.jl script to the matching .md, so that the prose,
# the code and the printed results can never drift apart: the code in the
# markdown is the code that produced the output beneath it.
#
#   julia --project=examples examples/generate.jl
#
# Script format: cells separated by a line of `#-`. Within a cell, lines
# starting with `#md` are prose and everything else is code that is shown and
# run. Whatever the code prints follows it as an output block.
#

const SCRIPTS = ["01-why-sectors", "02-sector-schemes", "03-category-rules", "04-multipass", "05-kriging"]

function capture(f)
  original = stdout
  rd, wr = redirect_stdout()
  reader = @async read(rd, String)
  try
    f()
  finally
    redirect_stdout(original)
    close(wr)
  end
  fetch(reader)
end

isprose(line) = startswith(line, "#md")
unprose(line) = rstrip(length(line) > 4 ? line[5:end] : "")

function render(name)
  # the working tree may be CRLF; strip CR so the cell separator regex matches
  source = filter(!isequal('\r'), read(joinpath(@__DIR__, name * ".jl"), String))
  sandbox = Module(Symbol(name))
  # shared plumbing is preloaded rather than shown: every script may use
  # drillholes() and the little reporting helpers without repeating them
  Base.include(sandbox, joinpath(@__DIR__, "drillholes.jl"))
  Base.include(sandbox, joinpath(@__DIR__, "report.jl"))
  out = IOBuffer()

  for cell in split(source, r"^#-+$"m)
    lines = split(strip(cell, '\n'), '\n')
    prose = [unprose(l) for l in lines if isprose(l)]
    code = strip(join([l for l in lines if !isprose(l)], '\n'), '\n')

    isempty(prose) || println(out, join(prose, '\n'), '\n')
    isempty(code) && continue

    println(out, "```julia\n", code, "\n```\n")
    printed = capture(() -> include_string(sandbox, code, name * ".jl"))
    isempty(strip(printed)) || println(out, "```\n", rstrip(printed), "\n```\n")
  end

  target = joinpath(@__DIR__, name * ".md")
  write(target, rstrip(String(take!(out))) * '\n')
  target
end

for name in SCRIPTS
  @info "rendering $name"
  println(stderr, "  wrote ", render(name))
end
