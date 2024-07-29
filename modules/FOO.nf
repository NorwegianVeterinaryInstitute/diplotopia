process FOO {
        debug true
        input:
        tuple val(thing)

        script:
        """
        cat $thing > text
        """
}



process FOO_PATH {
        input:
        path(necat_input)

        output:
        path(necat_input), emit: necat_input_ch

        script:
        """
        echo $necat_input
        """
}


