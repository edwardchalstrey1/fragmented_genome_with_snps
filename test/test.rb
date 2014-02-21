puts 'hello world'

require "rinruby"

myr = RinRuby.new(echo=false)
x = myr.pull "10"
myr.quit

puts x