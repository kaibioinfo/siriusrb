require "colorized_string"
module Sirius
  module Utils
    module Depiction

      PALETTE = [Java::JavaAwt::Color.new(166,206,227),Java::JavaAwt::Color.new(31,120,180),Java::JavaAwt::Color.new(178,223,138),Java::JavaAwt::Color.new(51,160,44),Java::JavaAwt::Color.new(251,154,153),Java::JavaAwt::Color.new(227,26,28),Java::JavaAwt::Color.new(253,191,111),Java::JavaAwt::Color.new(255,127,0),Java::JavaAwt::Color.new(202,178,214),Java::JavaAwt::Color.new(106,61,154),Java::JavaAwt::Color.new(255,255,153),Java::JavaAwt::Color.new(177,89,40)].to_a.each_slice(2).to_a.select {|x| x.size==2}.transpose.flatten
      COLOR_NAMES = ColorizedString.colors.drop(2).each_slice(2).to_a.select {|x| x.size==2}.transpose.flatten
      class DepictionPanel < Java::JavaxSwing::JPanel
        attr_accessor :img
        def initialize(img)
          super()
          @img = img
          setPreferredSize(Java::JavaAwt::Dimension.new(img.width, img.height))
        end
        def paint(g)
          g.clearRect(0, 0, width, height);
          g.drawImage(@img, 0,0,nil);
        end
      end
      class DepictionWindow < Java::JavaxSwing::JFrame
        def initialize(bufferedImage, w, h)
          super()
          imgPanel = DepictionPanel.new(bufferedImage)
          imgPanel.setPreferredSize(Java::JavaAwt::Dimension.new(w, h))
          scrollbar = Java::JavaxSwing::JScrollPane.new(imgPanel)
          getContentPane.setLayout(Java::JavaAwt::BorderLayout.new())
          getContentPane.add(scrollbar, Java::JavaAwt::BorderLayout::CENTER)
          setDefaultCloseOperation(Java::JavaxSwing::WindowConstants::DISPOSE_ON_CLOSE)
          setPreferredSize(Java::JavaAwt::Dimension.new(w+160,h))
          pack()
          setVisible(true)
        end
      end
      class DepictionWindowTable < Java::JavaxSwing::JFrame
        def initialize(images)
          super()
          imgPanel = DepictionPanel.new(images.first)
          list = Java::JavaxSwing::JList.new(images.to_java)
          list.setCellRenderer do |alist, avalue,anindex,issel,hasfoc|
            imgPanel.img = avalue
            imgPanel.setBackground(anindex.even? ? Java::JavaAwt::Color::WHITE : Java::JavaAwt::Color::LIGHT_GRAY )
            imgPanel
          end
          list.setPrototypeCellValue(images.first)
          scrollbar = Java::JavaxSwing::JScrollPane.new(list)
          getContentPane.setLayout(Java::JavaAwt::BorderLayout.new())
          getContentPane.add(scrollbar, Java::JavaAwt::BorderLayout::CENTER)
          setDefaultCloseOperation(Java::JavaxSwing::WindowConstants::DISPOSE_ON_CLOSE)
          pack()
          setVisible(true)
        end
      end

      Options2Colors = {
        :violet => Java::JavaAwt::Color.new(142, 68, 173),
        :red => Java::JavaAwt::Color.new(231, 76, 60 ),
        :blue => Java::JavaAwt::Color.new(52, 152, 219),
        :orange => Java::JavaAwt::Color.new(243, 156, 18 ),
        :green => Java::JavaAwt::Color.new( 39, 174, 96 )
      }

      extend self
      include_package "org.openscience.cdk.depict"
      # def quick_depict_many molecules, options
      #   # render at max 4 molecules next to each other
      #   size = java.awt.Toolkit.getDefaultToolkit().getScreenSize()
      #   minWidthPerMolecule = 620
      #   margin = 160
      #   maxcols = ((size.getWidth-margin)/minWidthPerMolecule).floor        
      #   maxrows = (molecules.size/maxcols.to_f).ceil
      #   windowSize = [maxcols*minWidthPerMolecule, maxrows*400]
      #   if maxrows == 1 && maxcols > molecules.size
      #     windowSize = [minWidthPerMolecule * molecules.size, 400]
      #   end
      #   ############################################
      #   gen = DepictionGenerator.new()
      #   Options2Colors.map {|key, value|
      #     if options[key]
      #       gen = gen.withHighlight(options[key], value)
      #     end
      #   }
      #   if options[:glow]
      #     gen = gen.withOuterGlowHighlight
      #   end
      #   depiction = gen.#withSize(windowSize[0], windowSize[1]).
      #   #withFillToFit.
      #   withZoom(1.5).
      #   withAromaticDisplay.withAtomColors.depict(molecules, maxrows, maxcols)
      #   image = depiction.toImg
      #   DepictionWindow.new(image, windowSize[0], windowSize[1])
      # end

      def quick_depict_many molecules, options
        if options[:table_layout]
          return quick_depict_table(molecules,options)
        end
        # render at max 4 molecules next to each other
        size = java.awt.Toolkit.getDefaultToolkit().getScreenSize()
        maxcols = [molecules.size, 4].min      
        maxrows = (molecules.size/maxcols.to_f).ceil
        ############################################
        gen = DepictionGenerator.new()
        gen = gen.withMolTitle
        Options2Colors.map {|key, value|
          if options[key]
            gen = gen.withHighlight(options[key], value)
          end
        }
        if options[:glow]
          gen = gen.withOuterGlowHighlight
        end
        depiction = gen.withZoom(1.5).withAromaticDisplay.withAtomColors.depict(molecules, maxrows, maxcols)
        image = depiction.toImg
        DepictionWindow.new(image, image.getWidth, image.getHeight)
      end
      def quick_depict_table moltab, options
        cols = moltab.map(&:size).max
        rows = moltab.size
        gen = DepictionGenerator.new()
        gen = gen.withMolTitle.withBackgroundColor(Java::JavaAwt::Color.new(0,0,0,0))
        d=Java::JavaAwt::Toolkit.getDefaultToolkit().getScreenSize();
        w=(d.getWidth-240)/cols
        h=(d.getHeight-120)/cols
        w=options[:width] if options[:width]
        h=options[:height] if options[:height]
        gen = gen.withZoom(1.5).withAromaticDisplay.withAtomColors.withSize(w,h)
        rows=moltab.map {|t|
          Options2Colors.map {|key, value|
          if options[key]
            gen = gen.withHighlight(options[key], value)
          end
          }
          if options[:glow]
            gen = gen.withOuterGlowHighlight
          end
          gen.depict(t, 1, cols).toImg
        }
        DepictionWindowTable.new(rows)
      end

    end
  end
end