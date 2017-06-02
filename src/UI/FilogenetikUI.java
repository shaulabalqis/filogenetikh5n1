/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package UI;
import UI.HasilFilogenetik;
import filogenetik.Filogenetik;
import globalalignment.globalalignment;
/**
 *
 * @author Kiky
 */
public class FilogenetikUI extends javax.swing.JFrame {

    /**
     * Creates new form FilogenetikUI
     */
    public FilogenetikUI() {
	initComponents();
	eSequence.setVisible(false);
	eJumlah.setVisible(false);
	eGap.setVisible(false);
	eMatriks.setVisible(false);
    }

    /**
     * This method is called from within the constructor to initialize the form.
     * WARNING: Do NOT modify this code. The content of this method is always
     * regenerated by the Form Editor.
     */
    @SuppressWarnings("unchecked")
    // <editor-fold defaultstate="collapsed" desc="Generated Code">//GEN-BEGIN:initComponents
    private void initComponents() {

        korea = new javax.swing.JCheckBox();
        thailand = new javax.swing.JCheckBox();
        china = new javax.swing.JCheckBox();
        france = new javax.swing.JCheckBox();
        germany = new javax.swing.JCheckBox();
        jumlah = new javax.swing.JTextField();
        teksJumlah = new javax.swing.JLabel();
        teksSequence = new javax.swing.JLabel();
        proses = new javax.swing.JButton();
        batal = new javax.swing.JButton();
        teksGap = new javax.swing.JLabel();
        gap = new javax.swing.JTextField();
        teksMatriks = new javax.swing.JLabel();
        PAM250 = new javax.swing.JCheckBox();
        BLOSUM62 = new javax.swing.JCheckBox();
        eJumlah = new javax.swing.JLabel();
        eSequence = new javax.swing.JLabel();
        eGap = new javax.swing.JLabel();
        eMatriks = new javax.swing.JLabel();

        setDefaultCloseOperation(javax.swing.WindowConstants.EXIT_ON_CLOSE);
        setTitle("Filogenetik UPGMA");

        korea.setText("Korea");
        korea.setEnabled(false);
        korea.addMouseListener(new java.awt.event.MouseAdapter() {
            public void mouseClicked(java.awt.event.MouseEvent evt) {
                koreaMouseClicked(evt);
            }
        });

        thailand.setText("Thailand");
        thailand.setEnabled(false);
        thailand.addMouseListener(new java.awt.event.MouseAdapter() {
            public void mouseClicked(java.awt.event.MouseEvent evt) {
                thailandMouseClicked(evt);
            }
        });

        china.setText("China");
        china.setEnabled(false);
        china.addMouseListener(new java.awt.event.MouseAdapter() {
            public void mouseClicked(java.awt.event.MouseEvent evt) {
                chinaMouseClicked(evt);
            }
        });

        france.setText("France");
        france.setEnabled(false);
        france.addMouseListener(new java.awt.event.MouseAdapter() {
            public void mouseClicked(java.awt.event.MouseEvent evt) {
                franceMouseClicked(evt);
            }
        });

        germany.setText("Germany");
        germany.setEnabled(false);
        germany.addMouseListener(new java.awt.event.MouseAdapter() {
            public void mouseClicked(java.awt.event.MouseEvent evt) {
                germanyMouseClicked(evt);
            }
        });

        jumlah.setText("jumlah");
        jumlah.addMouseListener(new java.awt.event.MouseAdapter() {
            public void mouseClicked(java.awt.event.MouseEvent evt) {
                jumlahMouseClicked(evt);
            }
        });
        jumlah.addFocusListener(new java.awt.event.FocusAdapter() {
            public void focusLost(java.awt.event.FocusEvent evt) {
                jumlahFocusLost(evt);
            }
        });
        jumlah.addKeyListener(new java.awt.event.KeyAdapter() {
            public void keyReleased(java.awt.event.KeyEvent evt) {
                jumlahKeyReleased(evt);
            }
        });

        teksJumlah.setText("masukkan jumlah sequence yang akan diproses (maksimal 5, minimal 3)");

        teksSequence.setText("pilih sequence yang akan diproses");

        proses.setText("proses");
        proses.addMouseListener(new java.awt.event.MouseAdapter() {
            public void mouseClicked(java.awt.event.MouseEvent evt) {
                prosesMouseClicked(evt);
            }
        });

        batal.setText("batal");
        batal.addMouseListener(new java.awt.event.MouseAdapter() {
            public void mouseClicked(java.awt.event.MouseEvent evt) {
                batalMouseClicked(evt);
            }
        });

        teksGap.setText("masukkan gap penjajaran");
        teksGap.setToolTipText("");

        gap.setText("gap");
        gap.addMouseListener(new java.awt.event.MouseAdapter() {
            public void mouseClicked(java.awt.event.MouseEvent evt) {
                gapMouseClicked(evt);
            }
        });
        gap.addFocusListener(new java.awt.event.FocusAdapter() {
            public void focusLost(java.awt.event.FocusEvent evt) {
                gapFocusLost(evt);
            }
        });

        teksMatriks.setText("pilih jenis matriks untuk penjajaran");

        PAM250.setText("PAM250");
        PAM250.addMouseListener(new java.awt.event.MouseAdapter() {
            public void mouseClicked(java.awt.event.MouseEvent evt) {
                PAM250MouseClicked(evt);
            }
        });

        BLOSUM62.setText("BLOSUM62");
        BLOSUM62.addMouseListener(new java.awt.event.MouseAdapter() {
            public void mouseClicked(java.awt.event.MouseEvent evt) {
                BLOSUM62MouseClicked(evt);
            }
        });

        eJumlah.setForeground(new java.awt.Color(255, 0, 0));
        eJumlah.setText("jumlah sequence tidak boleh kosong");

        eSequence.setForeground(new java.awt.Color(255, 0, 0));
        eSequence.setText("sequence harus dipilih sebanyak jumlah sequence");
        eSequence.setToolTipText("");

        eGap.setForeground(new java.awt.Color(255, 0, 0));
        eGap.setText("gap tidak boleh kosong");

        eMatriks.setForeground(new java.awt.Color(255, 0, 0));
        eMatriks.setText("matriks belum dipilih");

        javax.swing.GroupLayout layout = new javax.swing.GroupLayout(getContentPane());
        getContentPane().setLayout(layout);
        layout.setHorizontalGroup(
            layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(javax.swing.GroupLayout.Alignment.TRAILING, layout.createSequentialGroup()
                .addContainerGap(javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                .addComponent(batal)
                .addGap(18, 18, 18)
                .addComponent(proses)
                .addGap(80, 80, 80))
            .addGroup(layout.createSequentialGroup()
                .addGap(32, 32, 32)
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addComponent(teksJumlah)
                    .addGroup(layout.createSequentialGroup()
                        .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                            .addComponent(france)
                            .addComponent(thailand)
                            .addComponent(teksSequence)
                            .addComponent(korea)
                            .addComponent(china)
                            .addComponent(germany, javax.swing.GroupLayout.PREFERRED_SIZE, 163, javax.swing.GroupLayout.PREFERRED_SIZE)
                            .addComponent(eSequence))
                        .addGap(65, 65, 65)
                        .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                            .addGroup(layout.createSequentialGroup()
                                .addGap(1, 1, 1)
                                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                                    .addGroup(layout.createSequentialGroup()
                                        .addComponent(gap, javax.swing.GroupLayout.PREFERRED_SIZE, 34, javax.swing.GroupLayout.PREFERRED_SIZE)
                                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                                        .addComponent(eGap))
                                    .addComponent(teksGap)))
                            .addComponent(PAM250)
                            .addComponent(BLOSUM62)
                            .addComponent(teksMatriks)
                            .addComponent(eMatriks)))
                    .addGroup(layout.createSequentialGroup()
                        .addComponent(jumlah, javax.swing.GroupLayout.PREFERRED_SIZE, 69, javax.swing.GroupLayout.PREFERRED_SIZE)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(eJumlah)))
                .addContainerGap(122, Short.MAX_VALUE))
        );
        layout.setVerticalGroup(
            layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(layout.createSequentialGroup()
                .addGap(28, 28, 28)
                .addComponent(teksJumlah)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(jumlah, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(eJumlah))
                .addGap(18, 18, 18)
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING, false)
                    .addGroup(layout.createSequentialGroup()
                        .addComponent(teksGap)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                            .addComponent(gap, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                            .addComponent(eGap))
                        .addGap(31, 31, 31)
                        .addComponent(teksMatriks)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(PAM250)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                        .addComponent(BLOSUM62)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(eMatriks))
                    .addGroup(layout.createSequentialGroup()
                        .addComponent(teksSequence)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(korea)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(thailand)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(china)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(france)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(germany)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(eSequence, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED, 38, Short.MAX_VALUE)
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(proses)
                    .addComponent(batal))
                .addGap(31, 31, 31))
        );

        pack();
    }// </editor-fold>//GEN-END:initComponents

    private void jumlahMouseClicked(java.awt.event.MouseEvent evt) {//GEN-FIRST:event_jumlahMouseClicked
        eJumlah.setVisible(false);
	if (jumlah.getText().equals("jumlah")){
	    jumlah.setText(null);
	 }
    }//GEN-LAST:event_jumlahMouseClicked

    private void jumlahFocusLost(java.awt.event.FocusEvent evt) {//GEN-FIRST:event_jumlahFocusLost
       if (jumlah.getText().equals("")){
	   jumlah.setText("jumlah");
       }
    }//GEN-LAST:event_jumlahFocusLost

    private void jumlahKeyReleased(java.awt.event.KeyEvent evt) {//GEN-FIRST:event_jumlahKeyReleased
         if (jumlah.getText().equals("5") || jumlah.getText().equals("4")|| jumlah.getText().equals("3")){
	    korea.setEnabled(true);
	    thailand.setEnabled(true);
	    china.setEnabled(true);
	    france.setEnabled(true);
	    germany.setEnabled(true);
	} 
	else {
	    korea.setSelected(false);
	    thailand.setSelected(false);
	    china.setSelected(false);
	    france.setSelected(false);
	    germany.setSelected(false);
	    korea.setEnabled(false);
	    thailand.setEnabled(false);
	    china.setEnabled(false);
	    france.setEnabled(false);
	    germany.setEnabled(false);
	}
    }//GEN-LAST:event_jumlahKeyReleased

    private void gapMouseClicked(java.awt.event.MouseEvent evt) {//GEN-FIRST:event_gapMouseClicked
        eGap.setVisible(false);
	if (gap.getText().equals("gap")){
	    gap.setText(null);
	 }
    }//GEN-LAST:event_gapMouseClicked

    private void gapFocusLost(java.awt.event.FocusEvent evt) {//GEN-FIRST:event_gapFocusLost
        if (gap.getText().equals("")){
	   gap.setText("gap");
       }
    }//GEN-LAST:event_gapFocusLost

    private void PAM250MouseClicked(java.awt.event.MouseEvent evt) {//GEN-FIRST:event_PAM250MouseClicked
        eMatriks.setVisible(false);
	PAM250.setSelected(true);
	BLOSUM62.setSelected(false);
    }//GEN-LAST:event_PAM250MouseClicked

    private void BLOSUM62MouseClicked(java.awt.event.MouseEvent evt) {//GEN-FIRST:event_BLOSUM62MouseClicked
        eMatriks.setVisible(false);
	BLOSUM62.setSelected(true);
	PAM250.setSelected(false);
    }//GEN-LAST:event_BLOSUM62MouseClicked

    private void prosesMouseClicked(java.awt.event.MouseEvent evt) {//GEN-FIRST:event_prosesMouseClicked
      if (jumlah.getText().equals("")||jumlah.getText().equals("jumlah")){
	  eJumlah.setVisible(true);
      }
      else{
	  filogenetik.Filogenetik.jumlah=Integer.parseInt(jumlah.getText());
	  filogenetik.Filogenetik.init();
      }
      if (korea.isSelected()==false && thailand.isSelected()==false && china.isSelected()==false && france.isSelected()==false && germany.isSelected()==false){
	  eSequence.setVisible(true);
      }
      if (gap.getText().equals("")||gap.getText().equals("gap")){
	  eGap.setVisible(true);
      }
      else{
	  globalalignment.fungsi.gapInisialisasi=Integer.parseInt(gap.getText());
	  globalalignment.fungsi.gapIterasi=globalalignment.fungsi.gapInisialisasi;
      }
      if (PAM250.isSelected()==false && BLOSUM62.isSelected()==false){
	  eMatriks.setVisible(true);
      }
      else if (PAM250.isSelected()==true)
	  globalalignment.fungsi.pilihanmatriks=1;
      else
	  globalalignment.fungsi.pilihanmatriks=2;
      int[] temp=new int[5];
      int k=0;
      if (korea.isSelected()==true){
	  temp[k]=0;
	  k++;
      }
      if (thailand.isSelected()==true){
	  temp[k]=1;
	  k++;
      }
      if (china.isSelected()==true){
	  temp[k]=2;
	  k++;
      }
      if (france.isSelected()==true){
	  temp[k]=3;
	  k++;
      }
      if (germany.isSelected()==true){
	  temp[k]=4;
	  k++;
      }
      if (k!=filogenetik.Filogenetik.jumlah)
	  eSequence.setVisible(true);
      else{
	  for (int i=0; i<filogenetik.Filogenetik.jumlah; i++){
	    filogenetik.Filogenetik.idFile[i]=filogenetik.Filogenetik.idFile0[temp[i]];
	    filogenetik.Filogenetik.kluster[i]=filogenetik.Filogenetik.kluster0[temp[i]];
	    filogenetik.Filogenetik.namaFile[i]=filogenetik.Filogenetik.namaFile0[temp[i]];	   
	  }
	filogenetik.Filogenetik.bacaSequence();
	filogenetik.Filogenetik.alignment();
	filogenetik.Filogenetik.cetakSistem();
	for (int i=1; i<filogenetik.Filogenetik.jumlah; i++)
		filogenetik.Filogenetik.klusteringUPGMA();
	//UI.HasilFilogenetik prim2 =new UI.HasilFilogenetik();
	//prim2.setVisible(true);
	this.dispose();
      }
	  
    }//GEN-LAST:event_prosesMouseClicked

    private void koreaMouseClicked(java.awt.event.MouseEvent evt) {//GEN-FIRST:event_koreaMouseClicked
        eSequence.setVisible(false);
    }//GEN-LAST:event_koreaMouseClicked

    private void thailandMouseClicked(java.awt.event.MouseEvent evt) {//GEN-FIRST:event_thailandMouseClicked
        eSequence.setVisible(false);
    }//GEN-LAST:event_thailandMouseClicked

    private void chinaMouseClicked(java.awt.event.MouseEvent evt) {//GEN-FIRST:event_chinaMouseClicked
        eSequence.setVisible(false);
    }//GEN-LAST:event_chinaMouseClicked

    private void franceMouseClicked(java.awt.event.MouseEvent evt) {//GEN-FIRST:event_franceMouseClicked
        eSequence.setVisible(false);
    }//GEN-LAST:event_franceMouseClicked

    private void germanyMouseClicked(java.awt.event.MouseEvent evt) {//GEN-FIRST:event_germanyMouseClicked
        eSequence.setVisible(false);
    }//GEN-LAST:event_germanyMouseClicked

    private void batalMouseClicked(java.awt.event.MouseEvent evt) {//GEN-FIRST:event_batalMouseClicked
	this.dispose();
    }//GEN-LAST:event_batalMouseClicked

    /**
     * @param args the command line arguments
     */
    public static void main(String args[]) {
	/* Set the Nimbus look and feel */
	//<editor-fold defaultstate="collapsed" desc=" Look and feel setting code (optional) ">
        /* If Nimbus (introduced in Java SE 6) is not available, stay with the default look and feel.
	 * For details see http://download.oracle.com/javase/tutorial/uiswing/lookandfeel/plaf.html 
	 */
	try {
	    for (javax.swing.UIManager.LookAndFeelInfo info : javax.swing.UIManager.getInstalledLookAndFeels()) {
		if ("Nimbus".equals(info.getName())) {
		    javax.swing.UIManager.setLookAndFeel(info.getClassName());
		    break;
		}
	    }
	} catch (ClassNotFoundException ex) {
	    java.util.logging.Logger.getLogger(FilogenetikUI.class.getName()).log(java.util.logging.Level.SEVERE, null, ex);
	} catch (InstantiationException ex) {
	    java.util.logging.Logger.getLogger(FilogenetikUI.class.getName()).log(java.util.logging.Level.SEVERE, null, ex);
	} catch (IllegalAccessException ex) {
	    java.util.logging.Logger.getLogger(FilogenetikUI.class.getName()).log(java.util.logging.Level.SEVERE, null, ex);
	} catch (javax.swing.UnsupportedLookAndFeelException ex) {
	    java.util.logging.Logger.getLogger(FilogenetikUI.class.getName()).log(java.util.logging.Level.SEVERE, null, ex);
	}
	//</editor-fold>

	/* Create and display the form */
	java.awt.EventQueue.invokeLater(new Runnable() {
	    public void run() {
		new FilogenetikUI().setVisible(true);
	    }
	});
    }
    // Variables declaration - do not modify//GEN-BEGIN:variables
    private javax.swing.JCheckBox BLOSUM62;
    private javax.swing.JCheckBox PAM250;
    private javax.swing.JButton batal;
    private javax.swing.JCheckBox china;
    private javax.swing.JLabel eGap;
    private javax.swing.JLabel eJumlah;
    private javax.swing.JLabel eMatriks;
    private javax.swing.JLabel eSequence;
    private javax.swing.JCheckBox france;
    private javax.swing.JTextField gap;
    private javax.swing.JCheckBox germany;
    private javax.swing.JTextField jumlah;
    private javax.swing.JCheckBox korea;
    private javax.swing.JButton proses;
    private javax.swing.JLabel teksGap;
    private javax.swing.JLabel teksJumlah;
    private javax.swing.JLabel teksMatriks;
    private javax.swing.JLabel teksSequence;
    private javax.swing.JCheckBox thailand;
    // End of variables declaration//GEN-END:variables
}