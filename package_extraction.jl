using Pkg

# Get the current environment status (installed packages)
status_info = Pkg.status()

# Extract package names and versions in a readable format
package_list = ["$(pkg.name)@$(pkg.version)" for pkg in status_info]

# Write the list to a file
open("package_list.txt", "w") do file
    for pkg in package_list
        println(file, pkg)
    end
end

println("Package list has been written to package_list.txt")
