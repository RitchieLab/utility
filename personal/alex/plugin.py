#!/usr/bin/env python

#package loki

# plugin architecture from:
# http://martyalchin.com/2008/jan/10/simple-plugin-framework/


class Plugin(type):
	def __init__(cls, name, bases, attrs):
		print "Plugin.py::__init__(%s,%s,%s,%s)" % (cls,name,bases,attrs)
		if not hasattr(cls, 'plugins'):
			# This branch only executes when processing the mount point itself.
			# So, since this is a new plugin type, not an implementation, this
			# class shouldn't be registered as a plugin. Instead, it sets up a
			# list where plugins can be registered later.
			cls.plugins = []
		else:
			# This must be a plugin implementation, which should be registered.
			# Simply appending it to the list is all that's needed to keep
			# track of it later.
			cls.plugins.append(cls)
#Plugin
