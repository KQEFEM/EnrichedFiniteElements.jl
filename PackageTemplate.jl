using PkgTemplates

tpl = Template(;
	user = "KQEFEM",
	authors = "Kieran Quaine",
	license = "Creative Commons Attribution-NonCommercial 4.0 International License (CC BY-NC 4.0)",
	dir = ".",
	plugins = [
		Git(; manifest = true, ssh = true),
		GitHubActions(; x86 = true),
		Codecov(),
		Documenter{GitHubActions}(),

		TravisCI(),
	],
)

generate("EnrichedFiniteElements", tpl)
